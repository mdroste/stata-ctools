{smcl}
{* *! version 1.4.0 20261127}{...}
{viewerjumpto "Syntax" "cfuzzyjoin##syntax"}{...}
{viewerjumpto "Description" "cfuzzyjoin##description"}{...}
{viewerjumpto "Options for link" "cfuzzyjoin##options_link"}{...}
{viewerjumpto "Options for match" "cfuzzyjoin##options_match"}{...}
{viewerjumpto "Options for group" "cfuzzyjoin##options_group"}{...}
{viewerjumpto "Remarks" "cfuzzyjoin##remarks"}{...}
{viewerjumpto "Examples" "cfuzzyjoin##examples"}{...}
{viewerjumpto "Stored results" "cfuzzyjoin##results"}{...}
{viewerjumpto "Acknowledgments" "cfuzzyjoin##acknowledgments"}{...}
{title:Title}

{phang}
{bf:cfuzzyjoin} {hline 2} C-accelerated fuzzy joins: record linkage, similarity
joins of two datasets or columns, and grouping of similar strings


{marker syntax}{...}
{title:Syntax}

{pstd}
Record linkage: the best match of each observation (as {cmd:reclink})

{p 8 17 2}
{cmd:cfuzzyjoin link}
{varlist}
{cmd:using}
{it:filename}{cmd:,}
{opt idm:aster(varname)}
{opt idu:sing(varname)}
{opt g:en(newvarname)}
[{it:link_options}]

{pstd}
Similarity join of two datasets: every pair of similar strings (as {cmd:matchit})

{p 8 17 2}
{cmd:cfuzzyjoin match}
{it:idmaster} {it:txtmaster}
{cmd:using}
{it:filename}{cmd:,}
{opt idu:sing(varname)}
{opt txtu:sing(varname)}
[{it:match_options}]

{pstd}
Similarity of two columns of the data in memory (as {cmd:matchit} without {cmd:using})

{p 8 17 2}
{cmd:cfuzzyjoin match}
{it:varname1} {it:varname2}
[{cmd:,} {it:match_options}]

{pstd}
Groups of similar strings (as {cmd:strgroup})

{p 8 17 2}
[{cmd:by} {varlist}{cmd::}]
{cmd:cfuzzyjoin group}
{varname}
{ifin}{cmd:,}
{opt gen:erate(newvarname)}
{opt thresh:old(#)}
[{it:group_options}]

{synoptset 28 tabbed}{...}
{synopthdr:link_options}
{synoptline}
{synopt:{opt wm:atch(numlist)}}weights of agreement, one per variable; default 1 each{p_end}
{synopt:{opt wnom:atch(numlist)}}weights of disagreement; default {opt wmatch()}{p_end}
{synopt:{opt orb:lock(varlist|none)}}compare only records that agree on at least one of these variables; default all of {varlist} when it has four or more{p_end}
{synopt:{opt req:uired(varlist)}}compare only records that agree on all of these variables{p_end}
{synopt:{opt exa:ctstr(varlist)}}string variables compared as equal or not, without the bigram comparator{p_end}
{synopt:{opt exc:lude(filename)}}leave out the pairs, and their records, already matched in {it:filename}{p_end}
{synopt:{opt mins:core(#)}}a match must score above {it:#}; default 0.6{p_end}
{synopt:{opt minb:igram(#)}}bigram score at which strings count as agreeing; default 0.6{p_end}
{synopt:{opt uv:arlist(varlist)}}names of the matching variables in {it:filename}, if they differ{p_end}
{synopt:{opt upr:efix(string)}}prefix of the matching variables brought from {it:filename}; default {cmd:U}{p_end}
{synopt:{opt _m:erge(newvarname)}}name of the merge variable; default {cmd:_merge}{p_end}
{synopt:{opt deb:ug(filename)}}save the pairs found (idmaster, score, idusing) to {it:filename}{p_end}
{synopt:{opt verbose}}display detailed timing information{p_end}
{synopt:{opt thr:eads(#)}}maximum number of threads to use{p_end}
{synoptline}

{synoptset 28 tabbed}{...}
{synopthdr:match_options}
{synoptline}
{synopt:{opt sim:ilmethod(simfunc)}}how strings become grams; default {cmd:bigram}{p_end}
{synopt:{opt s:core(scorefunc)}}{cmd:jaccard} (the default), {cmd:simple} or {cmd:minsimple}{p_end}
{synopt:{opt w:eights(wgtfunc)}}weight grams by their frequency: {cmd:noweights} (the default), {cmd:simple}, {cmd:root} or {cmd:log}{p_end}
{synopt:{opt wgtf:ile(filename)}}take the gram frequencies from {it:filename}{p_end}
{synopt:{opt t:hreshold(#)}}keep pairs scoring at least {it:#}; default 0.5 (two datasets only){p_end}
{synopt:{opt stopw:ordsauto}}ignore grams found in more than {opt swthreshold()} of the observations{p_end}
{synopt:{opt swt:hreshold(#)}}frequency of stopwords per observation; default 0.2{p_end}
{synopt:{opt di:agnose}}report the most frequent grams and what the index can save{p_end}
{synopt:{opt g:enerate(newvarname)}}name of the score variable; default {cmd:similscore}{p_end}
{synopt:{opt over:ride}}replace the data in memory even when they have unsaved changes{p_end}
{synopt:{opt ti:me}}display the date and time at the start and the end{p_end}
{synopt:{opt f:lag(#)}}accepted for compatibility (progress reports){p_end}
{synopt:{opt keepm:ata}}accepted for compatibility; there are no Mata objects to keep{p_end}
{synopt:{opt verbose}}display detailed timing information{p_end}
{synopt:{opt thr:eads(#)}}maximum number of threads to use{p_end}
{synoptline}

{synoptset 28 tabbed}{...}
{synopthdr:group_options}
{synoptline}
{synopt:{opt norm:alize(shorter|longer|none)}}divide the edit distance by the length of the shorter string (the default), of the longer one, or by nothing{p_end}
{synopt:{opt first}}match only strings with the same first character{p_end}
{synopt:{opt noclean}}do not trim leading and trailing blanks before comparing{p_end}
{synopt:{opt force}}accepted for compatibility; {cmd:cfuzzyjoin group} has no observation limit{p_end}
{synopt:{opt verbose}}display detailed timing information{p_end}
{synopt:{opt thr:eads(#)}}maximum number of threads to use{p_end}
{synoptline}

{pstd}
{it:simfunc} is one of {cmd:bigram}, {cmd:ngram,} {it:n}, {cmd:ngram_circ,} {it:n},
{cmd:token}, {cmd:cotoken}, {cmd:scotoken}, {cmd:firstgram,} {it:n},
{cmd:tkngram,} {it:n}, {cmd:soundex}, {cmd:soundex_nara}, {cmd:soundex_fk},
{cmd:soundex_ext}, {cmd:nysiis_fk}, {cmd:token_soundex}, or
{cmd:tokenwrap,} {cmd:"}{it:f}{cmd:"} with {it:f} one of {cmd:soundex},
{cmd:soundex_nara}, {cmd:soundex_fk}, {cmd:soundex_ext} and {cmd:nysiis_fk}.
{cmd:bigram}, {cmd:ngram}, {cmd:ngram_circ}, {cmd:token}, {cmd:tkngram} and
{cmd:token_soundex} take a last argument {cmd:1} to count each gram once per
string (matchit's unit flag), as in {cmd:similmethod(bigram, 1)}.


{marker description}{...}
{title:Description}

{* ctools-credit:begin (generated by scripts/sync_attribution.py; edit config/attribution.json)}{...}
{pstd}
{cmd:cfuzzyjoin} reimplements {cmd:reclink} 1.7 (January 14, 2010) by Michael
Blasnik; {cmd:matchit} 1.5.2 (May 2020) by Julio D. Raffo; and {cmd:strgroup}
1.0.4 (December 17, 2021) by Julian Reif. See
{help cfuzzyjoin##acknowledgments:Acknowledgments} for details.
{* ctools-credit:end}{...}

{pstd}
{cmd:cfuzzyjoin} joins records whose keys are similar rather than equal:
names spelled differently, addresses typed in different ways. Its three
subcommands do the work of three community commands, with their syntax,
options, results and messages, in compiled code that runs on all processors.

{pstd}
{cmd:cfuzzyjoin link} finds, for each observation in memory (the master data),
the observations of {it:filename} (the using data) most similar to it on the
variables of {varlist}, as {cmd:reclink} does. String variables are compared
with reclink's bigram comparator, numeric variables (and strings in
{opt exactstr()}) as equal or not, and the agreement of each variable counts
with its weights; the score of a pair goes from 0 to 1. Records equal on every
variable match with score 1. The using records with the highest score, when it
exceeds {opt minscore()}, are the match; all of them when they tie. The result
keeps every master observation and adds the score, {opt idusing()}, the using
data's variables (its matching variables with the prefix {cmd:U}) and a
{cmd:_merge} variable (3 matched, 1 not).

{pstd}
{cmd:cfuzzyjoin match} with {cmd:using} compares the strings {it:txtmaster} of
the data in memory with {opt txtusing()} of {it:filename}, as {cmd:matchit}
does: each string becomes a vector of grams (bigrams by default, or tokens,
n-grams, phonetic keys), and the similarity of two strings is a Jaccard-type
score of their vectors, optionally weighted by the gram frequencies. Every pair
scoring at least {opt threshold()} is kept. The result replaces the data in
memory: one observation per pair with {it:idmaster}, {it:txtmaster},
{opt idusing()}, {opt txtusing()} and the score. Without {cmd:using},
{cmd:cfuzzyjoin match} adds the score of the two string variables of each
observation.

{pstd}
{cmd:cfuzzyjoin group} puts similar values of a string variable into groups,
as {cmd:strgroup} does: two strings match when their Levenshtein edit distance,
divided by the length of the shorter string (see {opt normalize()}), is at most
{opt threshold()}, and the groups are the sets of strings joined by chains of
matches. Groups are numbered from 1 in the order of their first observation.
{cmd:by} forms groups within by-groups and numbers them by-group after
by-group.


{marker options_link}{...}
{title:Options for link}

{phang}
{opt idmaster(varname)} is required and names the variable that identifies the
master observations. It may be numeric or string, and must have no duplicates
and no missing values.

{phang}
{opt idusing(varname)} is required and names the variable of {it:filename} that
identifies its observations, with the same conditions.

{phang}
{opt gen(newvarname)} is required and names the new variable that holds the
matching score ({cmd:float}, displayed {cmd:%5.4f}). It is missing for master
observations without a match.

{phang}
{opt wmatch(numlist)} gives the weight of agreement of each variable of
{varlist}, at least 1; typical values are 1 to 20. A variable whose agreement
says much about a true match, such as a name, gets a large weight; one with
many common values, such as a city, a small one.

{phang}
{opt wnomatch(numlist)} gives the weight of disagreement of each variable.
A small value means that disagreement is common even among true matches
(phone numbers change). The default is {opt wmatch()}.

{phang}
{opt orblock(varlist|none)} compares a master observation only with the using
observations that equal it on at least one of these variables (or-blocking).
This saves time; it is the default when {varlist} has four or more variables.
{cmd:orblock(none)} compares every pair. Variables of {opt required()} are left
out of the or-blocks. A master observation whose or-blocking strings are all
empty has no candidates.

{phang}
{opt required(varlist)} compares a master observation only with the using
observations that equal it on all of these variables (and-blocking). They also
count in the score.

{phang}
{opt exactstr(varlist)} names string variables that agree only when equal; the
bigram comparator is not used for them.

{phang}
{opt exclude(filename)} names a dataset of earlier matches, with the
{opt idmaster()} and {opt idusing()} variables. Master and using observations
whose identifiers appear in a pair with both identifiers present are left out
of the search, so that {cmd:cfuzzyjoin link} can be run again with other
settings on what is left. They stay in the result without a match.

{phang}
{opt minscore(#)} is the score a match must exceed; default 0.6.

{phang}
{opt minbigram(#)} is the bigram score at which two strings count as agreeing;
default 0.6. The bigram score of a variable becomes its weights of agreement
and disagreement with a sharp change at this value.

{phang}
{opt uvarlist(varlist)} names the matching variables of {it:filename} when they
differ from {varlist}, in the same order.

{phang}
{opt uprefix(string)} is the prefix of the using matching variables brought
into the result; default {cmd:U}. A using variable named like a master
variable of {varlist} takes the prefix; other using variables that the master
data already have are not brought in.

{phang}
{opt _merge(newvarname)} names the merge variable; default {cmd:_merge}.

{phang}
{opt debug(filename)} saves the pairs found, {opt idmaster()}, the score and
{opt idusing()}, sorted by {opt idmaster()}, as {cmd:reclink} does.

{phang}
{opt verbose} displays the time each step took; {opt threads(#)} sets the
maximum number of threads, from 1 to 256 (see {help ctools##threads:ctools}).


{marker options_match}{...}
{title:Options for match}

{phang}
{opt idusing(varname)} and {opt txtusing(varname)} are required with
{cmd:using}: the numeric identifier and the string variable of
{it:filename}. {it:idmaster} must be numeric and {it:txtmaster} a string.
The identifiers need not be unique.

{phang}
{opt similmethod(simfunc)} sets how strings become grams. {cmd:bigram} (the
default) takes every two consecutive characters; {cmd:ngram,} {it:n} every
{it:n}; {cmd:ngram_circ,} {it:n} also wraps the end of the string around to its
start, so that word order matters less; {cmd:token} splits at blanks;
{cmd:cotoken} takes pairs of consecutive tokens and {cmd:scotoken} tokens and
their pairs; {cmd:firstgram,} {it:n} the first {it:n} characters of each
token; {cmd:tkngram,} {it:n} tokens and n-grams. {cmd:soundex} and
{cmd:soundex_nara} take Stata's soundex codes of the string, {cmd:soundex_fk}
the full soundex key, {cmd:soundex_ext} the growing prefixes of that key and
{cmd:nysiis_fk} the full NYSIIS key; {cmd:token_soundex} and {cmd:tokenwrap}
apply them to each token. Strings too short for one gram are one gram. All
methods count bytes and are case sensitive.

{phang}
{opt score(scorefunc)} sets the score of two gram vectors with dot product N
and squared norms D1 and D2: {cmd:jaccard} N/sqrt(D1*D2) (the default),
{cmd:simple} 2*N/(D1+D2), {cmd:minsimple} N/min(D1,D2), at most 1.

{phang}
{opt weights(wgtfunc)} weights each gram by its frequency f over both string
variables: {cmd:simple} 1/f, {cmd:root} 1/sqrt(f), {cmd:log} 1/(log(f)+1).
Frequent grams then count less.

{phang}
{opt wgtfile(filename)} takes the frequencies from {it:filename}, with a
string variable {cmd:grams} and a numeric variable {cmd:freq}; grams the file
lacks have frequency 1. With {opt stopwordsauto} the frequencies of the data
are added to those of the file, as {cmd:matchit} does. It applies only with
{opt weights()}.

{phang}
{opt threshold(#)} keeps the pairs whose score is at least {it:#}; default 0.5.
Only pairs that share a gram are scored, even with {cmd:threshold(0)}.

{phang}
{opt stopwordsauto} treats grams more frequent than {opt swthreshold()} times
the number of observations of both datasets as stopwords: they are left out of
the vectors. This saves time on common grams at the cost of some matches.

{phang}
{opt diagnose} lists the most frequent grams of each dataset and the grams that
cost the index the most comparisons.

{phang}
{opt generate(newvarname)} names the score variable, {cmd:similscore} by
default. With {cmd:using}, variables of the result whose names clash are
renamed by appending a number, as {cmd:matchit} does.

{phang}
{opt override} lets {cmd:cfuzzyjoin match} replace data in memory that have
unsaved changes.

{phang}
{opt time}, {opt flag()} and {opt keepmata} are {cmd:matchit}'s options for
timing, progress reports and Mata objects; {opt time} displays the date and
time, the others have no effect.


{marker options_group}{...}
{title:Options for group}

{phang}
{opt generate(newvarname)} is required and names the new {cmd:long} variable
of group numbers. Missing strings, and observations excluded by {it:if} or
{it:in}, get missing values.

{phang}
{opt threshold(#)} is required and greater than 0: the largest normalized edit
distance at which two strings match. As in {cmd:strgroup}, it is used rounded
to seven decimal places. With {cmd:normalize(none)} it is a number of edits.

{phang}
{opt normalize(shorter|longer|none)} divides the Levenshtein distance by the
length of the shorter string of the pair (the default), of the longer one, or
not at all. Lengths count bytes.

{phang}
{opt first} matches only strings that start with the same letter (in either
case); strings that start with something else form one more set. Groups are
numbered set by set: other characters first, then a to z.

{phang}
{opt noclean} compares the strings as they are. By default leading and
trailing blanks are removed first, and strings that are blank become missing.

{phang}
{opt force} is accepted for compatibility; {cmd:strgroup} needs it above
10,000 observations, {cmd:cfuzzyjoin group} does not.


{marker remarks}{...}
{title:Remarks}

{pstd}
{bf:Speed.} {cmd:link} scores the candidate pairs of each master observation
in parallel, from the distinct value combinations of the using variables and
with indexes for the blocking variables; {cmd:reclink} reads the using data
once per master observation. {cmd:match} builds an inverted index of the using
grams and accumulates the dot products of each master string in one pass over
the strings that share a gram with it. {cmd:group} compares each distinct
string only with strings whose length leaves room for a match, skips pairs
already in one group and pairs that a count of their characters rules out,
and computes the remaining distances with a bit-parallel algorithm; it needs
memory in proportion to the number of strings, where {cmd:strgroup} needs
it in proportion to its square.

{pstd}
{bf:Results.} {cmd:link} computes reclink's scores with reclink's rounding of
each step to {cmd:float}, and {cmd:match} matchit's scores: without weights
they are equal to the last bit; with weights they differ by at most 1e-13
(relative) in the validation suite, because sums run in another order. {cmd:group} gives strgroup's group numbers.

{pstd}
{bf:Order of observations.} {cmd:link} sorts the result by {opt idusing()}, as
{cmd:reclink} leaves it, with the master observations without a match last.
{cmd:match} lists the pairs by master observation and, within one, in the
order of the using observations; {cmd:matchit} leaves pairs within a master
observation in an arbitrary order. The variable types of {cmd:match}'s result
are those {help compress} leaves, as with {cmd:matchit}.

{pstd}
{bf:Differences from the original commands.}

{phang2}
{cmd:link}: when a using value combination that several using records share
matches several master observations, every master observation gets every
such record; {cmd:reclink}'s many-to-many merge pairs them off in an arbitrary
order instead. Master values are compared as values, where {cmd:reclink}
writes them into macros: master strings with double quotes work
({cmd:reclink} stops with r(132)), and equal numbers that are not integers
agree ({cmd:reclink} compares a {cmd:float} 0.1 with the text of its value
and scores the pair as disagreeing). Master strings still reach the bigram
comparator as {cmd:reclink}'s {cmd:syntax} passes them: trimmed, with runs of
blanks and tabs made one blank. {opt orblock()} with all or-blocking strings
empty gives no candidates where {cmd:reclink} stops with a syntax error.
{opt wnomatch()} weights must be numbers.

{phang2}
{cmd:match}: {it:filename} may omit the {cmd:.dta} extension. Gram frequencies
missing from {opt wgtfile()} count as 1 also with {cmd:using}, where
{cmd:matchit} stops with a Mata error. {opt diagnose} lists the grams that both
datasets have. Only the built-in similarity, score and weight functions are
available; {cmd:matchit} also runs Mata functions a user adds. When
{opt txtusing()} is named like the score variable, the score is named
{cmd:similscore1}; {cmd:matchit} never returns.

{phang2}
{cmd:group}: blanks are trimmed unless {opt noclean} is given, as
{cmd:strgroup}'s help file says; {cmd:strgroup} 1.0.4 never trims them (its
{opt noclean} has no effect). There is no 10,000-observation limit and no
warning about duplicate values. When no observation is selected, the new
variable is all missing ({cmd:strgroup} stops with r(198)). {cmd:strgroup}'s
plugin reads only the first 244 bytes of each string; {cmd:cfuzzyjoin group}
compares whole strings.

{pstd}
The data are limited to 2,147,483,647 observations; strings are compared up
to 2,045 bytes, and longer {cmd:strL} values on their first 2,045 bytes.


{marker examples}{...}
{title:Examples}

{pstd}Link firms by name and address, blocking on the ZIP code or the city:{p_end}
{phang2}{cmd:. use firms_a, clear}{p_end}
{phang2}{cmd:. cfuzzyjoin link name address zip city using firms_b, idmaster(id) idusing(rid) gen(score) wmatch(10 8 2 2) wnomatch(8 4 1 1)}{p_end}
{phang2}{cmd:. list name Uname score if _merge == 3 & score < .8}{p_end}

{pstd}A second pass without blocking, on what the first pass left:{p_end}
{phang2}{cmd:. keep if _merge == 3}{p_end}
{phang2}{cmd:. save pass1}{p_end}
{phang2}{cmd:. use firms_a, clear}{p_end}
{phang2}{cmd:. cfuzzyjoin link name address using firms_b, idmaster(id) idusing(rid) gen(score) orblock(none) exclude(pass1)}{p_end}

{pstd}Every pair of names with a bigram similarity of at least 0.7:{p_end}
{phang2}{cmd:. use names_a, clear}{p_end}
{phang2}{cmd:. cfuzzyjoin match id name using names_b, idusing(rid) txtusing(rname) threshold(0.7)}{p_end}

{pstd}Token similarity weighted by rarity, ignoring frequent words:{p_end}
{phang2}{cmd:. cfuzzyjoin match id name using names_b, idusing(rid) txtusing(rname) similmethod(token) weights(log) stopwordsauto}{p_end}

{pstd}The similarity of two columns:{p_end}
{phang2}{cmd:. cfuzzyjoin match name1 name2, similmethod(ngram, 3) generate(sim3)}{p_end}

{pstd}Group spellings of firm names within each state:{p_end}
{phang2}{cmd:. bysort state: cfuzzyjoin group name, generate(firm) threshold(0.2)}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:cfuzzyjoin link} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:r(N)}}number of observations in the result{p_end}
{synopt:{cmd:r(matched)}}master observations with a match{p_end}
{synopt:{cmd:r(perfect)}}master observations with a perfect match{p_end}
{synopt:{cmd:r(unmatched)}}master observations without a match{p_end}

{pstd}
{cmd:cfuzzyjoin match} with {cmd:using} stores {cmd:r(N)}, the number of
pairs. {cmd:cfuzzyjoin group} stores:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:r(N)}}number of strings grouped{p_end}
{synopt:{cmd:r(groups)}}number of groups{p_end}
{synopt:{cmd:r(distinct)}}number of distinct strings{p_end}


{title:Author}

{pstd}
Michael Droste{break}
{browse "https://github.com/mdroste/stata-ctools":github.com/mdroste/stata-ctools}


{* ctools-attribution:begin (generated by scripts/sync_attribution.py; edit config/attribution.json)}{...}
{marker acknowledgments}{...}
{title:Acknowledgments}

{pstd}
{cmd:cfuzzyjoin} reimplements three community commands: {cmd:reclink}, written
by Michael Blasnik, as {cmd:cfuzzyjoin link}; {cmd:matchit}, written by Julio
D. Raffo (World Intellectual Property Organization), as
{cmd:cfuzzyjoin match}; and {cmd:strgroup}, written by Julian Reif (University
of Illinois), as {cmd:cfuzzyjoin group}. The syntax, options, scores, results,
and messages of each subcommand follow its original.

{pstd}
{cmd:cfuzzyjoin} is pinned to these versions:

{p2colset 9 27 29 2}{...}
{p2col:Command}{cmd:reclink}{p_end}
{p2col:Authors}Michael Blasnik{p_end}
{p2col:Version}1.7{p_end}
{p2col:Date}January 14, 2010{p_end}
{p2col:Header line}{cmd:*! v 1.7 14-Jan-2010 M. Blasnik:  Record Linkage}{p_end}
{p2col:Source}SSC ({cmd:ssc install reclink}), distributed January 18,
2010.{p_end}
{p2col:Latest version?}Yes, as of November 25, 2026.{p_end}
{p2colreset}{...}

{p2colset 9 27 29 2}{...}
{p2col:Command}{cmd:matchit}{p_end}
{p2col:Authors}Julio D. Raffo{p_end}
{p2col:Version}1.5.2{p_end}
{p2col:Date}May 2020{p_end}
{p2col:Header line}{cmd:*! 1.5.2 J.D. Raffo May 2020}{p_end}
{p2col:Source}SSC ({cmd:ssc install matchit}), distributed May 20,
2020.{p_end}
{p2col:Latest version?}Yes, as of November 25, 2026.{p_end}
{p2colreset}{...}

{p2colset 9 27 29 2}{...}
{p2col:Command}{cmd:strgroup}{p_end}
{p2col:Authors}Julian Reif{p_end}
{p2col:Version}1.0.4{p_end}
{p2col:Date}December 17, 2021{p_end}
{p2col:Header line}{cmd:*! strgroup 1.0.4 17dec2021 by Julian Reif}{p_end}
{p2col:Source}SSC ({cmd:ssc install strgroup}), distributed August 22,
2023.{p_end}
{p2col:Latest version?}Yes, as of November 25, 2026.{p_end}
{p2colreset}{...}

{pstd}
Julian Reif's GitHub repository (github.com/reifjulian/strgroup) has a newer
{cmd:strgroup}, 1.0.5 (January 2026), with a faster plugin; SSC distributes
1.0.4, which {cmd:cfuzzyjoin group} follows.

{pstd}
{cmd:cfuzzyjoin} changes only when ctools is updated to newer versions of
these commands and tested against them. {cmd:which} shows the version of each
command installed on your machine. {help ctools##acknowledgments:help ctools}
lists the pinned versions of all ctools commands.

{pstd}
Please report bugs and questions about {cmd:cfuzzyjoin} at
{browse "https://github.com/mdroste/stata-ctools/issues":github.com/mdroste/stata-ctools/issues}.
The authors of {cmd:reclink}, {cmd:matchit}, and {cmd:strgroup} do not
maintain {cmd:cfuzzyjoin}. Undocumented differences between {cmd:cfuzzyjoin}
and the original commands are bugs in {cmd:cfuzzyjoin}.

{pstd}
{cmd:matchit} computes gram frequencies with {cmd:freqindex} 1.3.1 (April
2019), also by Julio D. Raffo and distributed on SSC
({cmd:ssc install freqindex}); {cmd:cfuzzyjoin match} follows its similarity
functions too.

{pstd}
{cmd:strgroup}'s plugin computes the Levenshtein distance with code adapted
from David Necas's python-Levenshtein module. {cmd:cfuzzyjoin group} computes
the same distance with the bit-parallel algorithm of Myers (1999), in the form
Hyyrö (2003) gives for the distance between whole strings, and a banded
dynamic program (Ukkonen 1985) beyond 64 bytes.

{pstd}
If you use {cmd:cfuzzyjoin} in published work, please cite the original
command:

{phang2}Blasnik, M. 2007. RECLINK: Stata module to probabilistically match
records. Statistical Software Components S456876, Boston College Department of
Economics. Revised January 18, 2010.{p_end}

{phang2}Raffo, J. 2015. MATCHIT: Stata module to match two datasets based on
similar text patterns. Statistical Software Components S457992, Boston College
Department of Economics. Revised May 20, 2020.{p_end}

{phang2}Raffo, J., and S. Lhuillery. 2009. How to play the "Names Game":
Patent retrieval comparing different heuristics. {it:Research Policy} 38(10):
1617-1627.{p_end}

{phang2}Reif, J. 2010. STRGROUP: Stata module to match strings based on their
Levenshtein edit distance. Statistical Software Components S457151, Boston
College Department of Economics. Revised August 22, 2023.{p_end}

{phang2}Myers, G. 1999. A fast bit-vector algorithm for approximate string
matching based on dynamic programming. {it:Journal of the ACM} 46(3):
395-415.{p_end}

{phang2}Hyyrö, H. 2003. A bit-vector algorithm for computing Levenshtein and
Damerau edit distances. {it:Nordic Journal of Computing} 10(1): 29-39.{p_end}

{phang2}Ukkonen, E. 1985. Algorithms for approximate string matching.
{it:Information and Control} 64(1-3): 100-118.{p_end}
{* ctools-attribution:end}{...}
