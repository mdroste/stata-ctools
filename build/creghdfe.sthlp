{smcl}
{* *! version 1.2.6 20261111}{...}
{viewerjumpto "Syntax" "creghdfe##syntax"}{...}
{viewerjumpto "Description" "creghdfe##description"}{...}
{viewerjumpto "Options" "creghdfe##options"}{...}
{viewerjumpto "Prediction" "creghdfe##prediction"}{...}
{viewerjumpto "Examples" "creghdfe##examples"}{...}
{viewerjumpto "Stored results" "creghdfe##results"}{...}
{title:Title}

{phang}
{bf:creghdfe} {hline 2} C-accelerated high-dimensional fixed effects regression


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:creghdfe}
{depvar}
[{indepvars}]
{ifin}
{weight}
[{cmd:,} {opt a:bsorb(varlist)} {it:options}]

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Model}
{synopt:{opt a:bsorb(absvars[, suboptions])}}categorical variables representing fixed effects to absorb, with optional heterogeneous slopes and {it:newvar}{cmd:=} targets; optional; suboptions {opt save:fe} (or {opt g:enerate}), {opt replace}, {opt noi:sily}{p_end}
{synopt:{opt noab:sorb}}absorb only the constant (the default without {opt absorb()}){p_end}
{synopt:{opt nocon:stant}}do not report {cmd:_cons}{p_end}
{synopt:{opt keepsin:gletons}}keep singleton groups{p_end}
{synopt:{opt dof:adjustments(doftype)}}degrees of freedom adjustment method{p_end}
{synopt:{opt group:var(newvar)}}save mobility group identifier{p_end}
{synopt:{opt g:roup(groupvar)}}one regression row per group (a patent), whose variables are constant within it{p_end}
{synopt:{opt i:ndividual(indvar)}}with {opt group()}: the individuals (inventors) of each group; {it:indvar} in {opt absorb()} absorbs individual fixed effects{p_end}
{synopt:{opt ag:gregation(mean|sum)}}how the individual effects add up within a group; default {cmd:mean}{p_end}

{syntab:SE/Robust}
{synopt:{opt vce(vcetype)}}variance estimator; {opt un:adjusted} (default; also {opt ols} or {opt conventional}), {opt r:obust}, {opt cl:uster} {it:clustvars}, {opt dk:raay} [{it:#}], or {opt boot:strap}[{cmd:,} {it:bootopts}]{p_end}
{synopt:{opt cl:uster(clustvars)}}synonym for {opt vce(cluster} {it:clustvars}{opt )}{p_end}
{synopt:{opt fevce(bootstrap[, feopts])}}bootstrap standard errors of the {opt savefe} fixed effects{p_end}
{synopt:{opt wild:bootstrap}[{opt (wildopts)}]}wild (cluster) bootstrap p-values and confidence intervals of the coefficients{p_end}

{syntab:Convergence}
{synopt:{opt tol:erance(#)}}convergence tolerance; default is {cmd:1e-8}{p_end}
{synopt:{opt iterate(#)}}maximum iterations; default is {cmd:16000}{p_end}
{synopt:{opt thr:eads(#)}}maximum number of threads to use{p_end}
{synopt:{opt nostand:ardize}}do not standardize variables before iteration{p_end}
{synopt:{opt tra:nsform(str)}}{cmd:sym}metric_kaczmarz (default), {cmd:kac}zmarz or {cmd:cim}mino, as in {cmd:reghdfe}{p_end}
{synopt:{opt accel:eration(str)}}{cmd:cg} (default), {cmd:sd}, {cmd:a}itken or {cmd:none}, as in {cmd:reghdfe}{p_end}
{synopt:{opt tech:nique(str)}, {opt precond:itioner(str)}}accepted as in {cmd:reghdfe}; creghdfe solves with MAP{p_end}
{synopt:{opt prune}, {opt pool:size(#)}, {opt compact}, {opt par:allel(#)}, {opt fastreg:ress}}accepted as in {cmd:reghdfe}; no effect{p_end}

{syntab:Residuals}
{synopt:{opt res:iduals}}save the residuals in {cmd:_reghdfe_resid}, replacing it{p_end}
{synopt:{opt res:iduals(newvar)}}save the residuals in {it:newvar}; also {opt resid2(newvar)}{p_end}

{syntab:Advanced}
{synopt:{opt quad}}use quad-precision accumulation in CG solver{p_end}

{syntab:Data recovery}
{synopt:{opt snap:shot}}use {help preserve} to restore data after an error or interruption{p_end}
{synopt:{opt fast}}skip data recovery after an error or interruption{p_end}

{syntab:Reporting}
{synopt:{opt verbose}}display progress information and timing; also {opt verbose(#)} with {it:#} > 0{p_end}
{synopt:{opt l:evel(#)}}confidence level; and the other {help estimation options:display options}{p_end}
{synopt:{opt nohe:ader}}suppress the header{p_end}
{synopt:{opt notab:le}}suppress the coefficient table{p_end}
{synopt:{opt nofoot:note}}suppress the absorbed degrees of freedom table{p_end}
{synopt:{opt nowild:table}}suppress the wild bootstrap table{p_end}
{synopt:{opt nosamp:le}}do not save {cmd:e(sample)}{p_end}
{synopt:{opt time:it}}show {cmd:reghdfe}'s timer tables; {cmd:r(timer)}{p_end}
{synopt:{opt nowarn}}omit the warning of {opt keepsingletons}{p_end}
{synoptline}

{pstd}
{it:weight}s are allowed; {opt aweight}s, {opt fweight}s, and {opt pweight}s are supported.
{opt pweight}s imply {opt vce(robust)} unless clustering is specified.


{marker description}{...}
{title:Description}

{pstd}
{cmd:creghdfe} estimates linear regressions with multiple categorical fixed
effects. The fixed effects are removed by iterative projection in C. Without
{opt absorb()}, the command estimates a linear regression with a constant.
Factor-variable and time-series notation are allowed for the regressors, and,
as in {cmd:reghdfe}, for the dependent variable when it names one variable
({cmd:L.y} or the indicator {cmd:1.foreign}; {cmd:i.foreign} is r(9)).

{pstd}
Observations with missing model variables are excluded. Singleton fixed-effect
groups are dropped repeatedly until none remain (unless {opt keepsingletons}).
{cmd:e(sample)} marks the observations retained for estimation. Without
regressors, the model is the fixed effects and the constant.

{pstd}
Collinear regressors are omitted as {cmd:reghdfe} omits them: a regressor
whose partialled sum of squares is at most min(1e-6, {opt tolerance()}/10)
times its uncentered sum of squares is collinear with the fixed effects (with
a note). Among the remaining regressors, collinear variables are identified
using Mata's {cmd:invsym()} applied to the cross-product matrix in the
specified variable order. A dependent variable that passes the first rule is
perfectly explained by the fixed effects: as in {cmd:reghdfe} it is set to
zero, with a warning, so the coefficients and {cmd:e(rss)} are zero and
{cmd:e(F)} missing. A coefficient with zero estimate and zero variance is
marked omitted.

{pstd}
A saturated model (no residual degrees of freedom) is estimated as in
{cmd:reghdfe}: {cmd:e(df_r)} is 0, {cmd:e(V)} is zero, and {cmd:e(F)} and
{cmd:e(r2_a)} are missing. The F statistic is the Wald test of all retained
slopes under the chosen VCE, missing when the covariance matrix is rank
deficient.

{pstd}
Typing {cmd:creghdfe} without arguments replays the results with the same
options as estimation. {cmd:creghdfe} accepts {cmd:reghdfe}'s solver and memory
options (see {it:Convergence} below). It does not run {cmd:reghdfe}'s historical
versions: {opt version(3)} and {opt version(5)} are r(198) ({opt version()}
with any other number is r(9), as in {cmd:reghdfe}). Nor does it have
{cmd:reghdfe}'s programming options {opt nopartialout}, {opt noregress},
{opt keepmata} and {opt varlist_is_touse} (r(198)): they hand {cmd:reghdfe}'s
Mata object {cmd:HDFE} to other commands, and creghdfe, a compiled plugin, has
none.


{marker options}{...}
{title:Options}

{dlgtab:Model}

{phang}
{opt absorb(absvars[, suboptions])} specifies one or more categorical variables
whose fixed effects are to be absorbed. Each term is a numeric or string
variable, an {cmd:i.}{it:var} (or {cmd:ib}{it:#}{cmd:.}{it:var},
{cmd:ibn.}{it:var}: the base level does not matter), or an interaction
{it:var1}{cmd:#}{it:var2} absorbed as one effect. As in {cmd:reghdfe}, a term may also carry
heterogeneous slopes: {it:id}{cmd:#c.}{it:x} absorbs a slope on {it:x} for
each level of {it:id} (without an intercept), {it:id}{cmd:##c.}{it:x} an
intercept and a slope, and {it:id}{cmd:##c.(}{it:x z}{cmd:)} or
{it:id}{cmd:##(c.}{it:x} {cmd:c.}{it:z}{cmd:)} an intercept and several slopes.
When no absorbed effect has an intercept (slopes alone, as in
{cmd:absorb(}{it:id}{cmd:#c.}{it:x}{cmd:)}), nothing absorbs the constant: as
in {cmd:reghdfe}, the model is fitted without one, {cmd:_cons} is not
reported, {cmd:e(tss)} is uncentered and {cmd:e(report_constant)} is 0; such a
model needs a regressor. Each intercept and each
slope is one row of the absorbed degrees of freedom table and of
{cmd:e(extended_absvars)}.

{pmore}
A term may be preceded by {it:newvar}{cmd:=} to save its estimates in
{it:newvar}; its slopes are saved in {it:newvar}{cmd:Slope1},
{it:newvar}{cmd:Slope2}, ... The optional {opt savefe} suboption (or its
synonym {opt generate}) saves every effect not so named, in {cmd:__hdfe1__},
{cmd:__hdfe2__}, etc. (and {cmd:__hdfe1__Slope1}, ... for slopes), numbered by
term. A saved slope is the coefficient of the level, not its product with the
slope variable. As in {cmd:reghdfe}, every run first drops the variables
{cmd:__hdfe*} (so that a rerun, or {helpb bootstrap}, with {opt savefe} works);
the other output names must be new, except that suboption {opt replace}
replaces an existing {it:newvar} target ({it:newvar}{cmd:Slope}{it:#} must
still be new). Suboption {opt noisily} shows the absorbed effects as parsed;
any other suboption is r(198). Individual fixed effects depend on the normalization used.
They are calculated iteratively, subject to {opt iterate()}, and sum to y - xb -
residual; with slopes, each effect's intercept and slopes are fitted together
within each level, as {cmd:reghdfe} does. On a poorly connected graph (low
mobility between, say, workers and firms) the individual levels are poorly
determined in both commands. Named targets cannot be combined with
{opt fevce()}; use {opt savefe}.

{phang}
{opt dofadjustments(list)} selects how redundant fixed-effect coefficients are
counted, as in {cmd:reghdfe}. {it:list} may be {opt all} (default), {opt none},
or any of {opt first:pair} or {opt pair:wise}, {opt cl:usters} and
{opt cont:inuous}. {opt pairwise} counts each effect's redundant coefficients
from the mobility groups it forms with every earlier effect; {opt firstpair}
examines the first two effects only. {opt clusters} treats effects nested
within the cluster variable as costing no degrees of freedom and sets them
aside before mobility groups are counted. {opt all} is {opt pairwise clusters}
{opt continuous}. {opt continuous} counts a heterogeneous slope as redundant in
each level where its variable (divided by its standard deviation) is
constant, or zero when its effect has no intercept, to within 1e-6. Counts
that may be too low (every slope, and intercepts after the first pair) are
marked {cmd:?} in the table and in column {cmd:Inexact?} of
{cmd:e(dof_table)}.

{phang}
{opt group(groupvar)}, {opt individual(indvar)} and {opt aggregation(mean|sum)}
fit {cmd:reghdfe}'s models of group outcomes with individual fixed effects. The
data have one observation per group and individual (a patent and each of its
inventors); the dependent variable, the regressors, the cluster variables and
the weight must be constant within each group (r(498) otherwise), and a group
and an individual must identify an observation (r(459)). The regression has one
row per group, {cmd:e(N)} groups. With {opt individual()}, {it:indvar} must be
one of the absorbed effects; it enters each group's row as the mean of the
effects of the group's individuals, or with {cmd:aggregation(sum)} their sum.
Each individual costs one degree of freedom, none counted as redundant (marked
{cmd:?}), as in {cmd:reghdfe}. As there, a group with an individual found in no
other group is a singleton. The projection is solved with LSMR, as
{cmd:reghdfe}'s {cmd:technique(lsmr)}, 100 times more tightly than
{opt tolerance()}. Not available with {opt individual()}: fweights, heterogeneous
slopes, saving the effects ({opt savefe}, {it:newvar}{cmd:=}),
{opt vce(bootstrap)}, {opt fevce()}, {opt wildbootstrap()} and {opt groupvar()}.
With aweights or pweights the estimates are the weighted least-squares fit (the
regression on explicit individual-share variables), which {cmd:reghdfe} 6.13.1
does not reach. {opt group()} alone keeps one row per group and absorbs the other
effects as usual.

{phang}
{opt groupvar(newvar)} creates a new variable containing the mobility group
identifier for each observation. Mobility groups are connected components in the
network formed by the first two fixed effects that are not nested within the
cluster variable, as in {cmd:reghdfe}. As there, no variable is created when no
such pair exists or degrees of freedom are not adjusted for mobility groups
({cmd:dofadjustments(none)} or {cmd:clusters} alone). With one effect, all
retained observations receive 1.

{dlgtab:SE/Robust}

{phang}
{opt vce()} reads its type as {cmd:reghdfe} does (case matters): any word
starting with {cmd:r} is {opt robust}, with {cmd:un} or {cmd:ols}
{opt unadjusted} (as is {cmd:conventional}), with {cmd:cl} {opt cluster} and
with {cmd:dk} {opt dkraay}; any other type is r(9).

{phang}
{opt vce(robust)} computes heteroskedasticity-robust standard errors.

{phang}
{opt vce(cluster clustvars)} computes cluster-robust standard errors, clustering
on the specified variables, each of which may be a string variable, an
{cmd:i.}{it:var}, or an interaction {it:var1}{cmd:#}{it:var2}. With several
variables (up to ten), the standard errors are multi-way clustered as in
{cmd:reghdfe} (Cameron, Gelbach and Miller 2011): the clustered covariances of
every intersection of the cluster variables are added (an odd number of
variables) or subtracted (an even number), with one small-sample adjustment
using the fewest clusters, and negative eigenvalues of the result are set to
zero (with a warning) while the coefficient block keeps its own correction.
That combination can leave {cmd:e(V)} indefinite, and {cmd:ereturn post} may
then set it to zero, as it does for {cmd:reghdfe} with some orders of the
cluster variables. {cmd:creghdfe} then posts a positive semidefinite
{cmd:e(V)} with the same coefficient block instead, changing only the
constant's row and column.
{cmd:e(df_r)} is the fewest clusters minus one (at most the unclustered
{cmd:e(df_r)}). An absorbed effect nested in any cluster variable costs no
degrees of freedom. {opt wildbootstrap} takes one cluster variable.

{phang}
{opt vce(dkraay} [{it:#}]{opt )} computes Driscoll-Kraay standard errors,
robust to heteroskedasticity and to correlation across panels and over time,
as {cmd:reghdfe} does. The data must be {helpb tsset} or {helpb xtset} with a time
variable. The moments are summed within each time period and a Bartlett kernel
of bandwidth {it:#} (lags + 1; default floor(4(T/100)^(2/9)) + 1 for T periods)
weights their autocovariances; a bandwidth above T is reduced to T, with a
warning. The time variable is treated as the cluster variable for the degrees
of freedom: {cmd:e(df_r)} is at most T - 1, and absorbed effects nested in the
periods cost none. {cmd:e(dkraay_bw)} is the bandwidth used.

{phang}
With {cmd:pweight}s the standard errors are robust unless {opt vce(cluster)} is
specified; {opt vce(unadjusted)} and {opt vce(ols)} are an error, as in {cmd:reghdfe}.

{phang}
{opt vce(bootstrap} [{cmd:,} {it:bootopts}]{opt )} computes bootstrap standard
errors (a pairs bootstrap: observations or clusters are resampled with
replacement). {cmd:reghdfe} has no such option; the replicates run inside the
plugin on the data loaded for the fit, much faster than the {helpb bootstrap}
prefix. The results are those {cmd:vce(bootstrap)} posts after official
estimation commands, so {helpb bootstrap postestimation:estat bootstrap} reports
the normal, percentile and bias-corrected intervals. Inference is z-based;
{cmd:e(df_r)} is not posted, and the header reports a Wald chi-squared test of
the regressors. The fit's own (unadjusted) covariance is posted as
{cmd:e(V_modelbased)}. A replicate whose omitted regressors differ from the full
sample's is missing, as with the {cmd:bootstrap} prefix. {it:bootopts} are:

{p2colset 9 32 34 2}{...}
{p2col:{opt r:eps(#)}}replications; default 50{p_end}
{p2col:{opt seed(#)}}{cmd:set seed} {it:#} before the first draw{p_end}
{p2col:{opt cl:uster(clustvar)}}resample clusters (one variable, {cmd:i.}{it:var}, or {it:var1}{cmd:#}{it:var2}){p_end}
{p2col:{opt str:ata(varlist)}}resample within strata; a cluster is a unit within each stratum{p_end}
{p2col:{opt draws(stata|fast)}}how Stata's random numbers become draws; default {cmd:stata}{p_end}
{p2col:{opt mse}}center the variance on {cmd:e(b)} instead of the replicate mean ({cmd:e(vcetype)} is then {cmd:Bstrap *} and {cmd:e(mse)} is {cmd:mse}, as with {cmd:bootstrap}){p_end}
{p2col:{opt l:evel(#)}}level of the posted intervals and of the table at estimation; the command's own {opt level()} sets it too, and the two may not differ{p_end}
{p2col:{opt nodots}, {opt dots(#)}}progress display{p_end}
{p2col:{opt sav:ing(filename}[{cmd:,} {opt replace} {opt double} {opt every(#)}]{opt )}}save the replicates' coefficients as {helpb bootstrap} saves them{p_end}
{p2colreset}{...}

{pmore}
{opt saving()} writes one observation per replicate and one variable per term
of {cmd:e(b)}, named, labeled and with the characteristics {helpb bootstrap}
gives them ({cmd:_b_}{it:name}, or {cmd:_bs_}{it:#}), so that {helpb bstat}
reads the file. A replicate that failed is missing throughout; omitted terms
are zero. The values are {cmd:float} unless {opt double} is given (as with
{cmd:bootstrap}); {cmd:e(V)} is computed from the replicates in double either
way.

{pmore}
Resampling acts on the estimation sample, after singletons are dropped, as the
{cmd:bootstrap} prefix resamples {cmd:e(sample)}. Each replicate drops its own
singletons and judges collinearity again. Random numbers come from Stata's
generator: {cmd:set seed}, {opt seed()}, {cmd:set rng} and {cmd:set rngstream}
apply, {cmd:e(rngstate)} records the state before the first draw, and the
generator is left past the draws taken. With {cmd:draws(stata)}, the default,
the draws are exactly those of {helpb bsample}, so {cmd:creghdfe} reproduces
{cmd:bootstrap _b, reps(}{it:#}{cmd:) seed(}{it:#}{cmd:): reghdfe} {it:...} for the same
seed. Its uniforms (two per resampled observation or cluster and replicate)
are drawn in Stata and passed to the plugin in chunks of replicates.
{cmd:draws(fast)} takes a master seed from Stata's generator and draws in the
plugin: reproducible from the seed, faster for large observation-level
bootstraps, but not {cmd:bsample}'s draws. Weights are allowed (official {cmd:vce(bootstrap)}
refuses them): {cmd:aweight}s and {cmd:pweight}s multiply the resampled rows,
and without {opt cluster()} {cmd:fweight}s are resampled as the observations
they stand for. With {cmd:draws(stata)} these are {cmd:bsample}'s draws on the
fweight-expanded data, in the order {helpb expand} leaves it, so the results
are those of the {cmd:bootstrap} prefix after {cmd:expand}. With
{opt strata()} and observations resampled, {cmd:draws(stata)} also takes the
uniforms that {cmd:bsample}, as the prefix calls it, takes after its draws
(one per observation, without changing the sample), so that later replicates
and the generator's final state match the prefix too. Heterogeneous slopes in
{opt absorb()} are re-estimated in every replicate, as the fit estimates
them.

{phang}
{opt fevce(bootstrap} [{cmd:,} {it:feopts}]{opt )} computes bootstrap standard
errors of the fixed effects that {opt absorb(..., savefe)} saves, by the same
pairs bootstrap. It requires {opt savefe}. Each replicate re-estimates the model
on its resample and recovers the effects of the levels it contains; a level's
standard error is the standard deviation of its effect over the replicates in
which it is estimable. The standard error of the {it:g}th effect is saved in
{cmd:__fese}{it:g}{cmd:__} on the estimation sample (a slope's, saved by
{opt savefe} in {cmd:__hdfe}{it:g}{cmd:__Slope}{it:h}, in
{cmd:__fese}{it:g}{cmd:__Slope}{it:h}), missing where fewer than two
replicates estimated the level. A slope is estimable in a level where the
replicate's slope variable varies; where the full fit estimates a level's
slope but a replicate cannot, that replicate does not estimate the level's
intercept either. With {opt vce(bootstrap)} the coefficients
and the fixed effects share one set of replicates; resampling options given to
both must agree, and may be given to either, and {opt level()}, {opt mse} and
{opt dots()} are those of {opt vce()} ({opt fevce()} may repeat them). Alone, {opt fevce()} leaves
{cmd:e(V)} as {opt vce()} computes it. {it:feopts} are {it:bootopts}
({opt saving()} saves the coefficients' replicates, and with
{opt vce(bootstrap)} it goes there) and:

{p2colset 9 32 34 2}{...}
{p2col:{opt sestub(stub)}}name the standard errors {it:stub}{cmd:1}, {it:stub}{cmd:2}, ...{p_end}
{p2col:{opt nreps}}also save the number of replicates that estimated each level, in {cmd:__fenrep}{it:g}{cmd:__} (or {it:stub}{it:g}{cmd:_nrep}){p_end}
{p2col:{opt dse(newvar)}}save the standard error of the summed fixed effects d, the sum over the absorbed effects{p_end}
{p2col:{opt ci}}also save percentile bounds at {opt level()} in {cmd:__felb}{it:g}{cmd:__} and {cmd:__feub}{it:g}{cmd:__} (or {it:stub}{it:g}{cmd:_lb} and {cmd:_ub}), as {cmd:e(ci_percentile)} is computed; every replicate's effects are kept in memory (4 bytes each), and r(909) is issued when they do not fit{p_end}
{p2colreset}{...}

{pmore}
A level is estimable in a replicate when it was drawn and the connected set of
levels it belongs to in the full sample did not split; a replicate whose omitted
regressors differ from the full sample's estimates no fixed effects.
{cmd:e(fe_bs_table)} reports, per effect, the levels, how many were never and
always estimable, and the minimum and median count of replicates that estimated
a level.

{pmore}
With two or more absorbed effects, individual fixed effects are identified only
up to a normalization within each connected set, and their bootstrap standard
errors depend on that normalization; they can be large where the sets are
thinly connected. {cmd:creghdfe} then displays a warning (also on replay), and
{cmd:e(fe_norm_warning)} is 1. The summed effect d is identified, so its
standard errors from {opt dse()} do not depend on the normalization.

{phang}
{opt wildbootstrap}[{opt (wildopts)}] adds wild bootstrap p-values and
confidence intervals for the coefficients, after the analytic fit: e(b),
e(V) and the coefficient table are unchanged, and a second table reports,
for each tested coefficient, the bootstrap p-value of H0: b = 0 and the
interval found by inverting the bootstrap test. It needs {opt vce(robust)}
or {opt vce(cluster} {it:clustvar}{opt )}; with clusters it is the wild
cluster bootstrap, and with {opt vce(robust)} it draws one weight per
observation. By default the null hypothesis is imposed (WCR, the restricted
wild cluster bootstrap). The algorithm and its conventions are those of
{browse "https://ideas.repec.org/c/boc/bocode/s458121.html":boottest}
(Roodman, MacKinnon, Nielsen and Webb 2019): for the same {opt seed()} the
weights, and hence the p-values, are boottest's. {it:wildopts} are:

{p2colset 9 32 34 2}{...}
{p2col:{opt r:eps(#)}}replications; default 999{p_end}
{p2col:{opt seed(#)}}{cmd:set seed} {it:#} before the draws{p_end}
{p2col:{opt weight:type(name)}}{cmd:rademacher} (default), {cmd:mammen}, {cmd:webb}, {cmd:normal} or {cmd:gamma}{p_end}
{p2col:{opt nonull}}do not impose the null (WCU){p_end}
{p2col:{opt coef:s(names)}}coefficients to test; default every estimated regressor{p_end}
{p2col:{opt p:type(name)}}{cmd:symmetric} (default), {cmd:equaltail}, {cmd:lower} or {cmd:upper}{p_end}
{p2col:{opt noenum:erate}}draw weights even when all 2^C sign vectors fit in {opt reps()}{p_end}
{p2col:{opt fedgp(full|partialled)}}how absorbed effects not nested in the bootstrap clusters enter the bootstrap; default {cmd:full}{p_end}
{p2col:{opt noci}}skip the confidence intervals{p_end}
{p2col:{opt l:evel(#)}}level of the intervals{p_end}
{p2col:{opt bootcl:uster(varname)}}draw the weights at this (coarser) clustering instead{p_end}
{p2col:{opt nodots}, {opt dots(#)}}progress display{p_end}
{p2colreset}{...}

{pmore}
With Rademacher weights and 2^C <= {opt reps()} (C bootstrap clusters), all
2^C sign vectors are used instead of random draws, as in boottest: the
replications become 2^C and the random-number state is left unchanged. With
11 or fewer clusters a note suggests {cmd:weighttype(webb)}. Otherwise the
weights come from Stata's generator: {cmd:set seed}, {opt seed()},
{cmd:set rng} and {cmd:set rngstream} apply, and {cmd:e(wb_rngstate)}
records the state before the draws. Stata's {helpb wildbootstrap} draws
differently, so its results are not reproduced draw for draw.

{pmore}
Absorbed effects nested in the bootstrap clusters (as with
{cmd:absorb(firm) vce(cluster firm)}) cost nothing: each replication is a
few sums over clusters. An effect that crosses the clusters (year effects
with firm clusters) is partialled out of each replicate's perturbed
residuals: the bootstrap of the full model, as boottest after
{cmd:reghdfe} or {cmd:regress} with dummies computes it.
{cmd:fedgp(partialled)} instead bootstraps the partialled regression, as
Stata's {cmd:wildbootstrap areg} does, at the cost of the nested case.
{cmd:e(wb_fe_path)} records which ({cmd:nested}, {cmd:projected} or
{cmd:partialled}).

{pmore}
P-values count the replications whose statistic is strictly beyond the
observed one, over the replications with a defined statistic, as boottest
does. The intervals are the set of values b0 whose p-value exceeds
1 - {opt level()}/100, bounded exactly (boottest stops its search within
one replication's step of the level, so its bounds can differ slightly).
A bound that does not exist within 1000 standard errors is missing, and a
set that is not an interval is reported by its outer bounds; the table
notes both. {opt wildbootstrap()} may be combined with
{opt fevce(bootstrap)}, but not with {opt vce(bootstrap)}. With
{opt vce(robust)} and {cmd:fweight}s it needs {opt bootcluster()}.

{dlgtab:Convergence}

{phang}
{opt tolerance(#)} specifies a positive convergence tolerance for the iterative
demeaning algorithm. The default is {cmd:1e-8}.

{phang}
{opt iterate(#)} specifies a positive maximum number of iterations for the
demeaning algorithm (and for recovering {opt savefe} effects). The default is
{cmd:16000}, as in {cmd:reghdfe}.

{phang}
{opt threads(#)} sets the maximum number of threads, from 1 to 256.
The default follows the {help ctools##threads:ctools thread setting}.

{phang}
{opt nostandardize} specifies that variables should not be standardized
before the iterative algorithm. This may affect convergence speed.

{phang}
{opt transform(str)} and {opt acceleration(str)} select the variant of the
method of alternating projections, with {cmd:reghdfe}'s names, abbreviations
and errors: transforms {cmd:symmetric_kaczmarz} (default), {cmd:kaczmarz} and
{cmd:cimmino}; accelerations {cmd:conjugate_gradient} ({cmd:cg}, default),
{cmd:steepest_descent} ({cmd:sd}), {cmd:aitken} and {cmd:none} ({cmd:off}).
{cmd:kaczmarz} with {cmd:cg} runs, with {cmd:reghdfe}'s warning that it is
unlikely to converge; {cmd:rand_kaczmarz} is r(3498), as in {cmd:reghdfe}.
The estimates agree with the default's to within the tolerance.

{phang}
{opt technique(map|lsmr|lsqr)} and {opt preconditioner(str)} are checked as
{cmd:reghdfe} checks them ({cmd:technique(cg)} is r(3498), as there), but
creghdfe always solves with the method of alternating projections (MAP), to
the same {opt tolerance()} (with {opt individual()}, with LSMR). The
estimates agree to within the tolerance; {cmd:e(ic)} counts MAP iterations.
{opt prune}, {opt poolsize(#)}, {opt compact}, {opt fastregress} and
{opt parallel(#)} ({cmd:reghdfe}'s memory, speed and worker-process options)
are accepted and have no effect: creghdfe holds the data in C and runs in
threads ({opt threads()}).

{phang}
{cmd:reghdfe}'s undocumented solver properties may be set as there:
{opt abort(0)} keeps the projection reached after {opt iterate()} iterations
when it has not converged, with {cmd:e(converged)} = 0 (and, as
{cmd:reghdfe}'s MAP solver reports it, {cmd:e(ic)} = 0), instead of r(430).
{opt miniter(#)}, {opt min_ok(#)}, {opt accel_start(#)} and
{opt accel_freq(#)} are accepted and have no effect (in {cmd:reghdfe},
{opt miniter()} applies to LSMR and LSQR only, and creghdfe's Aitken
acceleration already uses {cmd:reghdfe}'s defaults). Any other
{it:key}{cmd:(}{it:value}{cmd:)} is r(198), as in {cmd:reghdfe}.

{dlgtab:Advanced}

{phang}
{opt quad} uses quad-precision (double-double) accumulation in the conjugate
gradient solver for improved numerical accuracy. This is slower but may help
with ill-conditioned problems.

{dlgtab:Residuals}

{phang}
{opt residuals} (abbreviated {opt res} or {opt resid}) stores the residuals
in {cmd:_reghdfe_resid}, replacing that variable when it exists, as in
{cmd:reghdfe}, so that reruns and {helpb bootstrap} work.

{phang}
{opt residuals(newvar)} (also {opt res(newvar)}, {opt resid(newvar)} and
{opt resid2(newvar)}) stores the residuals in a new variable. Only one of
the two forms may be given (r(198)).

{dlgtab:Data recovery}

{pstd}
By default, an error or interruption removes new output variables and clears
{cmd:e()}. Failure to remove fixed effects within the iteration limit returns
error 430.

{phang}
{opt snapshot} uses {help preserve} to copy the whole dataset for recovery.
Stata may write the copy to disk when it exceeds {cmd:c(max_preservemem)}.

{phang}
{opt fast} skips data recovery. An error or interruption can leave partially
filled output variables. It cannot be combined with {opt snapshot}.


{dlgtab:Reporting}

{phang}
{opt verbose} displays progress and elapsed times. {cmd:reghdfe}'s
{opt verbose(#)} is accepted: a positive level does the same.

{phang}
{opt level(#)} and the other {help estimation options:display options}
({opt noomitted}, {opt vsquish}, {opt baselevels}, {opt cformat()},
{opt coeflegend}, ...) format the coefficient table. {opt noheader},
{opt notable} and {opt nofootnote} suppress the header, the table and the table
of absorbed degrees of freedom. These options also apply on replay.

{phang}
{opt keepsingletons} retains singleton groups. Keeping these groups can bias
inference; a warning says so unless {opt nowarn} (or {opt verbose(-1)}) is
specified.

{phang}
{opt nosample} does not save {cmd:e(sample)}, as in {cmd:reghdfe}.

{phang}
{opt timeit} shows {cmd:reghdfe}'s two timer tables (top-level steps and
partialling out), filled with creghdfe's steps, and leaves the second in
{cmd:r(timer)}. Unlike {cmd:reghdfe}, it leaves Stata's {helpb timer}s alone.


{marker prediction}{...}
{title:Prediction}

{pstd}
{cmd:predict newvar} defaults to {opt xb}, the linear prediction excluding
absorbed effects. {opt stdp} gives its standard error. Estimate with
{opt residuals} or {opt residuals()} to use the following prediction options:

{p 8 12 2}{opt residuals}: stored regression residuals.{p_end}
{p 8 12 2}{opt d}: sum of absorbed effects, y - xb - residual.{p_end}
{p 8 12 2}{opt xbd}: fitted values including absorbed effects, y - residual.{p_end}
{p 8 12 2}{opt dresiduals}: absorbed effects plus residual, y - xb.{p_end}

{pstd}
Residual and fixed-effect output variables are missing outside the estimation
sample. Their names, including those from {opt groupvar()} and {opt savefe},
must be distinct new variables.

{marker examples}{...}
{title:Examples}

{pstd}Setup:{p_end}
{phang2}{cmd:. sysuse auto, clear}{p_end}

{pstd}Basic regression with a single fixed effect:{p_end}
{phang2}{cmd:. creghdfe price mpg weight, absorb(foreign)}{p_end}

{pstd}With robust standard errors:{p_end}
{phang2}{cmd:. creghdfe price mpg weight, absorb(foreign) vce(robust)}{p_end}

{pstd}Using panel data with two-way fixed effects:{p_end}
{phang2}{cmd:. webuse nlswork, clear}{p_end}
{phang2}{cmd:. creghdfe ln_wage age ttl_exp tenure, absorb(idcode year)}{p_end}

{pstd}With clustered standard errors:{p_end}
{phang2}{cmd:. creghdfe ln_wage age ttl_exp tenure, absorb(idcode year) vce(cluster idcode)}{p_end}

{pstd}Two-way clustering, and Driscoll-Kraay standard errors:{p_end}
{phang2}{cmd:. creghdfe ln_wage age ttl_exp tenure, absorb(idcode year) vce(cluster idcode year)}{p_end}
{phang2}{cmd:. xtset idcode year}{p_end}
{phang2}{cmd:. creghdfe ln_wage age ttl_exp tenure, absorb(idcode) vce(dkraay 3)}{p_end}

{pstd}Individual intercepts and experience slopes, saved as {cmd:fe} and {cmd:feSlope1}:{p_end}
{phang2}{cmd:. creghdfe ln_wage age tenure, absorb(fe=idcode##c.ttl_exp year)}{p_end}

{pstd}Save residuals:{p_end}
{phang2}{cmd:. creghdfe ln_wage age ttl_exp, absorb(idcode) resid2(myresid)}{p_end}

{pstd}With verbose output showing timing:{p_end}
{phang2}{cmd:. creghdfe ln_wage age ttl_exp tenure, absorb(idcode year) verbose}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:creghdfe} stores the following in {cmd:e()}:

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(df_m)}}model degrees of freedom{p_end}
{synopt:{cmd:e(df_r)}}residual degrees of freedom{p_end}
{synopt:{cmd:e(r2)}}R-squared{p_end}
{synopt:{cmd:e(r2_within)}}within R-squared{p_end}
{synopt:{cmd:e(r2_a)}}adjusted R-squared{p_end}
{synopt:{cmd:e(rmse)}}root mean squared error{p_end}
{synopt:{cmd:e(rss)}}residual sum of squares{p_end}
{synopt:{cmd:e(mss)}}model sum of squares{p_end}
{synopt:{cmd:e(tss)}}total sum of squares{p_end}
{synopt:{cmd:e(tss_within)}}within total sum of squares{p_end}
{synopt:{cmd:e(df_a)}}absorbed degrees of freedom (adjusted){p_end}
{synopt:{cmd:e(df_a_initial)}}levels of the absorbed effects, before redundant coefficients are removed{p_end}
{synopt:{cmd:e(df_a_redundant)}}redundant coefficients of the absorbed effects, nested ones included{p_end}
{synopt:{cmd:e(df_a_nested)}}degrees of freedom from nested FEs{p_end}
{synopt:{cmd:e(F)}}F-statistic{p_end}
{synopt:{cmd:e(r2_a_within)}}adjusted within R-squared{p_end}
{synopt:{cmd:e(ll)}}log likelihood{p_end}
{synopt:{cmd:e(ll_0)}}log likelihood of the model with the fixed effects only{p_end}
{synopt:{cmd:e(sumweights)}}sum of the weights ({cmd:e(N)} unweighted){p_end}
{synopt:{cmd:e(ic)}}iterations of the fixed-effect projection{p_end}
{synopt:{cmd:e(converged)}}1: the projection converged; 0 after {opt abort(0)} without convergence{p_end}
{synopt:{cmd:e(N_full)}}observations before singletons were dropped, {cmd:e(N)} + {cmd:e(num_singletons)}{p_end}
{synopt:{cmd:e(report_constant)}}1 if {cmd:_cons} is reported (0 with {opt noconstant}, or without an absorbed intercept){p_end}
{synopt:{cmd:e(rank)}}number of retained slope coefficients{p_end}
{synopt:{cmd:e(N_hdfe)}}number of absorbed FE groups{p_end}
{synopt:{cmd:e(N_hdfe_extended)}}number of absorbed intercepts and slopes, as in {cmd:reghdfe} (the mobility groups are the second row's redundant count in {cmd:e(dof_table)}){p_end}
{synopt:{cmd:e(num_singletons)}}number of singleton observations dropped{p_end}
{synopt:{cmd:e(drop_singletons)}}1 if singletons are dropped (0 with {opt keepsingletons}){p_end}
{synopt:{cmd:e(N_clust)}}fewest clusters of any cluster variable (if {cmd:vce(cluster)}; time periods if {cmd:vce(dkraay)}){p_end}
{synopt:{cmd:e(N_clust}{it:#}{cmd:)}}number of clusters of cluster variable {it:#} (if {cmd:vce(cluster)} or {cmd:vce(dkraay)}){p_end}
{synopt:{cmd:e(N_clustervars)}}number of cluster variables (if {cmd:vce(cluster)} or {cmd:vce(dkraay)}){p_end}
{synopt:{cmd:e(dkraay_bw)}}Driscoll-Kraay bandwidth (if {cmd:vce(dkraay)}){p_end}
{synopt:{cmd:e(fe_N_reps)}}replicates run for the fixed effects (if {cmd:fevce(bootstrap)}){p_end}
{synopt:{cmd:e(fe_split_reps)}}replicates in which a connected set of levels split (if {cmd:fevce(bootstrap)}){p_end}
{synopt:{cmd:e(fe_norm_warning)}}1 with two or more absorbed effects: the normalization caveat applies (if {cmd:fevce(bootstrap)}){p_end}
{synopt:{cmd:e(wb_reps)}, {cmd:e(wb_misreps)}}wild bootstrap replications, and those without a defined statistic (if {cmd:wildbootstrap}){p_end}
{synopt:{cmd:e(wb_enumerated)}}1 if all 2^C sign vectors were used (if {cmd:wildbootstrap}){p_end}
{synopt:{cmd:e(wb_N_clust)}, {cmd:e(wb_level)}}bootstrap clusters, and the level of the intervals (if {cmd:wildbootstrap}){p_end}

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:creghdfe}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(indepvars)}}regressors partialled out, as {cmd:reghdfe} lists them ({cmd:2bn.rep78}, ...; base levels left out), and {cmd:_cons} when it is reported{p_end}
{synopt:{cmd:e(absorb)}}absorbed fixed effect variables{p_end}
{synopt:{cmd:e(absvars)}}absorbed effects as specified ({cmd:_cons} without {opt absorb()}){p_end}
{synopt:{cmd:e(extended_absvars)}}absorbed intercepts and slopes, one per row of {cmd:e(dof_table)}: {it:id} and {it:id}{cmd:#c.}{it:x}, without {cmd:i.} prefixes{p_end}
{synopt:{cmd:e(title)}}{cmd:HDFE Linear regression}{p_end}
{synopt:{cmd:e(title2)}}number of absorbed effects{p_end}
{synopt:{cmd:e(dofmethod)}}degrees-of-freedom adjustments applied{p_end}
{synopt:{cmd:e(vce)}}VCE type: {cmd:ols} (unadjusted), {cmd:robust}, {cmd:cluster}, {cmd:dkraay}, or {cmd:bootstrap}{p_end}
{synopt:{cmd:e(vcetype)}}{cmd:Robust} (if {cmd:vce(robust)}, {cmd:vce(cluster)} or {cmd:vce(dkraay)}); {cmd:Bootstrap} with {cmd:vce(bootstrap)}, {cmd:Bstrap *} with its {opt mse}{p_end}
{synopt:{cmd:e(title3)}}{cmd:Statistics robust to heteroskedasticity} (if {cmd:vce(cluster)} or {cmd:vce(dkraay)}){p_end}
{synopt:{cmd:e(title4)}}{cmd:and autocorrelation (Driscoll-Kraay)} (if {cmd:vce(dkraay)}){p_end}
{synopt:{cmd:e(clustvar)}}names of the cluster variables (if {cmd:vce(cluster)}; the time variable if {cmd:vce(dkraay)}){p_end}
{synopt:{cmd:e(clustvar}{it:#}{cmd:)}}name of cluster variable {it:#} (if {cmd:vce(cluster)} or {cmd:vce(dkraay)}){p_end}
{synopt:{cmd:e(resid)}}name of residual variable (if {cmd:resid} or {cmd:resid2()} specified){p_end}
{synopt:{cmd:e(predict)}}program used to implement {cmd:predict}{p_end}
{synopt:{cmd:e(estat_cmd)}}{cmd:_creghdfe_estat}: {cmd:estat summarize}, {cmd:vce}, {cmd:ic} (with {cmd:e(df_m)} + 1 degrees of freedom, as {cmd:reghdfe}) and {cmd:bootstrap}{p_end}
{synopt:{cmd:e(footnote)}}{cmd:_creghdfe_footnote}: the absorbed degrees of freedom table, shown by prefix commands{p_end}
{synopt:{cmd:e(model)}}{cmd:ols}{p_end}
{synopt:{cmd:e(marginsnotok)}}predictions {cmd:margins} refuses: {cmd:Residuals SCore}{p_end}
{synopt:{cmd:e(wtype)}, {cmd:e(wexp)}}weight type and expression (if weighted){p_end}
{synopt:{cmd:e(fevce)}}{cmd:bootstrap} (if {cmd:fevce(bootstrap)}){p_end}
{synopt:{cmd:e(wb_method)}}{cmd:wcr} or {cmd:wcu} (if {cmd:wildbootstrap}){p_end}
{synopt:{cmd:e(wb_weighttype)}, {cmd:e(wb_ptype)}}weight type and p-value type{p_end}
{synopt:{cmd:e(wb_coefs)}, {cmd:e(wb_bootclust)}}tested coefficients, and the bootstrap cluster variable{p_end}
{synopt:{cmd:e(wb_fe_path)}}{cmd:nested}, {cmd:projected} or {cmd:partialled}{p_end}
{synopt:{cmd:e(wb_rngstate)}}random-number state before the wild bootstrap's draws{p_end}
{synopt:{cmd:e(fe_se)}, {cmd:e(fe_nrep)}, {cmd:e(fe_dse)}}variables holding the fixed effects' standard errors, replicate counts and d's standard errors (if {cmd:fevce(bootstrap)}){p_end}
{synopt:{cmd:e(fe_lb)}, {cmd:e(fe_ub)}}variables holding the percentile bounds (if {cmd:fevce(bootstrap, ci)}; scalar {cmd:e(fe_level)} their level){p_end}
{synopt:{cmd:e(fe_rngstate)}, {cmd:e(fe_draws)}}random-number state before the first draw, and {opt draws()} (if {cmd:fevce(bootstrap)}){p_end}
{synopt:{cmd:e(fe_cluster_bs)}, {cmd:e(fe_strata_bs)}}resampled clusters and strata (if {cmd:fevce(bootstrap)} with {opt cluster()} or {opt strata()}){p_end}

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix of the estimators{p_end}
{synopt:{cmd:e(V_modelbased)}}the fit's own covariance (if {cmd:vce(bootstrap)}){p_end}
{synopt:{cmd:e(b_bs)}, {cmd:e(reps)}, {cmd:e(bias)}, {cmd:e(se)}, {cmd:e(z0)}}bootstrap means, complete replicates, bias, standard errors and bias-correction constants (if {cmd:vce(bootstrap)}){p_end}
{synopt:{cmd:e(ci_normal)}, {cmd:e(ci_percentile)}, {cmd:e(ci_bc)}}bootstrap intervals (if {cmd:vce(bootstrap)}){p_end}
{synopt:{cmd:e(boot_status)}}replicates that were complete, collinear, empty, did not converge, or failed otherwise (if {cmd:vce(bootstrap)}){p_end}
{synopt:{cmd:e(wb_t)}, {cmd:e(wb_p)}}analytic t statistics and wild bootstrap p-values of the tested coefficients (if {cmd:wildbootstrap}){p_end}
{synopt:{cmd:e(wb_ci)}, {cmd:e(wb_cistatus)}}confidence intervals (missing bound: open), and per coefficient open lower, open upper, disjoint, empty, and failed (a trial value whose p-value could not be evaluated, treated as not rejected) (unless {cmd:noci}){p_end}
{synopt:{cmd:e(fe_bs_table)}}per effect: levels, never and always estimable, minimum and median replicates estimating a level (if {cmd:fevce(bootstrap)}){p_end}
{synopt:{cmd:e(dof_table)}}per absorbed intercept and slope: categories, redundant, coefficients, inexact and nested, named as in {cmd:reghdfe} (rows {cmd:1.}{it:id}, {cmd:1.}{it:a}{cmd:#1.}{it:b}, {cmd:1.}{it:id}{cmd:#c.}{it:x}; columns {cmd:Categories}, {cmd:Redundant}, {cmd:Num Coefs}, {cmd:Inexact?}, {cmd:Nested?}){p_end}

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}


{title:Author}

{pstd}
Michael Droste{break}
{browse "https://github.com/mdroste/stata-ctools":github.com/mdroste/stata-ctools}


{title:Acknowledgments}

{pstd}
{cmd:creghdfe} estimates fixed-effects regressions similar to Sergio Correia's
{cmd:reghdfe}. The original {cmd:reghdfe} package is available at
{browse "https://github.com/sergiocorreia/reghdfe"}.

{pstd}
See also: Correia, S. 2016. "Linear Models with High-Dimensional Fixed Effects:
An Efficient and Feasible Estimator." Working paper.
{browse "http://scorreia.com/research/hdfe.pdf"}


{title:Also see}

{psee}
Manual: {bf:[R] areg}

{psee}
Online: {browse "https://github.com/sergiocorreia/reghdfe":reghdfe} (if installed),
{help areg}, {help xtreg}, {help ctools}
{p_end}
