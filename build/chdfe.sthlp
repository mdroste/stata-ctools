{smcl}
{* *! version 1.2.3 20261030}{...}
{viewerjumpto "Syntax" "chdfe##syntax"}{...}
{viewerjumpto "Description" "chdfe##description"}{...}
{viewerjumpto "Options" "chdfe##options"}{...}
{viewerjumpto "Examples" "chdfe##examples"}{...}
{viewerjumpto "Stored results" "chdfe##results"}{...}
{title:Title}

{pstd}{bf:chdfe} {hline 2} Remove fixed effects from variables

{marker syntax}{...}
{title:Syntax}

{p 8 12 2}{cmd:chdfe} {varlist} {ifin} [{it:weight}],
{opt a:bsorb(absvars)} {opt gen:erate(prefix)} [{opt sample(newvar)} {it:options}]

{p 8 12 2}{cmd:chdfe} {varlist} {ifin} [{it:weight}],
{opt a:bsorb(absvars)} {opt clear} [{opt keepv:ars(varlist)} {opt keepid:s} {it:options}]

{pstd}{cmd:aweight}s, {cmd:fweight}s, and {cmd:pweight}s are allowed.
The variable list accepts factor variables and time-series operators.

{marker description}{...}
{title:Description}

{pstd}{cmd:chdfe} returns residuals after removing the effects listed in
{opt absorb()}, following the user-written {cmd:hdfe} command. It uses the C
projection routines shared with {help creghdfe}. It also reports the degrees
of freedom absorbed by the effects. It does not fit a regression.

{pstd}Rows with missing inputs, missing absorbed or cluster categories, or
missing or zero weights are excluded. Singleton categories are removed
repeatedly until no more remain. Generated residuals are missing outside
this sample. An empty sample returns error 2000.

{marker options}{...}
{title:Options}

{pstd}{opt absorb(absvars)} is required. It accepts categorical variables,
{cmd:i.var}, categorical interactions such as {cmd:firm#year}, and group-specific
slopes such as {cmd:firm#c.time}, {cmd:firm##c.time}, and {cmd:firm##c.(time age)}.
{cmd:#c.} removes slopes only; {cmd:##c.} also removes group intercepts.
Categorical keys may be strings; {cmd:strL} keys are limited to 2,045 bytes per
value. Fixed-effect coefficient output is not available.

{pstd}{opt generate(prefix)} creates one new double variable per input, named
by adding the prefix to its name. For expanded factor or time-series terms,
names use the {cmd:hdfe} convention: {cmd:__}, with dots replaced by {cmd:__}
and interactions by {cmd:_X_}. Base and omitted indicators are excluded.
Output names must be new. Original variables and observation order are preserved.

{pstd}{opt sample(newvar)} saves a zero-one sample indicator with {opt generate()}.

{pstd}{opt clear} replaces the dataset with the selected rows and transformed
variables. It also retains the source variables needed for expanded terms and
simple weight variables. {opt keepvars()} retains additional variables. {opt keepids} keeps
category IDs named {cmd:__ID1__}, {cmd:__ID2__}, and so on, and cluster IDs named
{cmd:__CL1__}, and so on. Specify exactly one of {opt generate()} and {opt clear}.
A failed call leaves the original dataset intact.

{pstd}{opt clustervars(varlist)} supplies cluster categories for the absorbed
degrees-of-freedom calculation. Categorical interactions are allowed. Effects
nested within a cluster do not count toward absorbed degrees of freedom when
cluster adjustments are enabled.

{pstd}{opt dofadjustments(list)} selects {cmd:pairwise}, {cmd:firstpair},
{cmd:clusters}, and {cmd:continuous} adjustments. The default, {cmd:all}, selects
{cmd:pairwise clusters continuous}. {cmd:none} disables these adjustments.
For more than two categorical effects, the result is generally an upper bound
on the degrees of freedom absorbed, as with {cmd:hdfe}.

{pstd}{opt groupvar(newvar)} saves connected-component IDs for the first pair of
nonnested categorical effects. It requires {cmd:pairwise} or {cmd:firstpair}
adjustments. No variable is created when there is no eligible pair.

{pstd}{opt keepsingletons} retains singleton categories.

{pstd}{opt tolerance(#)} sets convergence tolerance; the default is {cmd:1e-8}.
{opt maxiterations(#)} limits iterations; the default is 10,000 and {cmd:.}
removes the practical limit. Failure to converge returns error 430.

{pstd}{opt poolsize(#)} sets the number of variables projected together; the
default is 10. A smaller value reduces projection workspace.

{pstd}{opt transform(str)} selects {cmd:symmetric_kaczmarz} (default),
{cmd:kaczmarz}, or {cmd:cimmino}. Abbreviations are allowed.
{opt acceleration(str)} selects {cmd:conjugate_gradient} ({cmd:cg}, default),
{cmd:steepest_descent} ({cmd:sd}), {cmd:aitken}, or {cmd:none}.
Conjugate gradient requires a symmetric transform. Iteration paths and counts
may differ from {cmd:hdfe}.

{pstd}{opt threads(#)} sets the maximum worker count, from 1 to 256.
Zero or omission uses the {help ctools##threads:ctools setting}.

{pstd}{opt verbose}, {opt verbose(#)} with a positive number, and {opt timeit}
report load, projection, and write times. {cmd:verbose(0)} is quiet.
{cmd:chdfe, version} reports the ctools release and stores {cmd:e(version)}.

{marker examples}{...}
{title:Examples}

{phang2}{cmd:. chdfe sales capital, absorb(firm year) generate(within_)}{p_end}
{phang2}{cmd:. chdfe sales capital [aw=weight], absorb(firm year) generate(r_) sample(used)}{p_end}
{phang2}{cmd:. chdfe sales capital, absorb(firm##c.year) clear keepvars(firm year)}{p_end}

{marker results}{...}
{title:Stored results}

{pstd}{cmd:e(df_a)} holds absorbed degrees of freedom. {cmd:e(N_hdfe)} is the
number of absorb terms; {cmd:e(N_hdfe_extended)} counts intercepts and slopes
separately. {cmd:e(mobility)} is the total redundant coefficient count and
{cmd:e(M_due_to_nested)} is the count due to cluster nesting.

{pstd}For each expanded effect {it:j}, {cmd:e(Kj)} is its category count,
{cmd:e(Mj)} its redundant coefficient count, {cmd:e(Mj_exact)} and
{cmd:e(Mj_nested)} indicate exact and nested counts, and {cmd:e(Gj)} identifies
its original absorb term. {cmd:e(absvars)} and {cmd:e(extended_absvars)} contain
the original and expanded effects. {cmd:e(cmd)} is {cmd:chdfe}.

{pstd}Additional results are {cmd:e(N)}, {cmd:e(num_singletons)},
{cmd:e(iteration)}, {cmd:e(converged)}, and {cmd:e(dofmethod)}.

{title:Also see}
{pstd}{help creghdfe}, {help ctools}, {help hdfe}
