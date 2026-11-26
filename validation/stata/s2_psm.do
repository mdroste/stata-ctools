* s2_psm: named cases share this file, each do-call has its own local scope.
gettoken __case 0 : 0
macro shift
if "`__case'" == "" local __case "s2_psm"

if "`__case'" == "s2_psm" {
do validation/support/validate_setup.do
set seed 20261104

/* _s2psm_cmp runs psmatch2 and then cpsmatch with the same arguments (the
   command line follows the test name and a colon) on the same data and sort
   order. It compares the return code, every r() scalar and macro psmatch2
   returns, its r(table) (values, row and column names), and every variable
   psmatch2 creates (values, variable labels with cpsmatch: for psmatch2:,
   the storage type of _id and _n#), plus char _weight[Type]. Values agree
   within 1e-9 (reldif); both missing is agreement. The reference is
   psmatch2 4.1.0 (see docs/ATTRIBUTION.md), held in S2REF_<variable>,
   scalars S2R_<name> and matrix S2T_ref.

   One documented departure (docs/COMPATIBILITY.md, "cpsmatch matching"):
   for ties with neighbor(k) or ai(k), k > 1, with replacement,
   _s2psm_ties_rule replaces the affected part of that reference with the
   rule cpsmatch documents, so those calls are compared with the rule and
   not with psmatch2. */

/* psmatch2's ties adds the donors tied with the nearest one and then takes
   the k - 1 next positions around the nearest: a tied donor is listed (and
   weighted) twice, and donors tied at the k-th distance are left out.
   cpsmatch matches every donor at most as far as the k-th nearest, each
   once. The functions below compute that rule by brute force from
   psmatch2's _treated, _pscore and _support and store what follows from it
   in the reference: _nn, _weight and _<outcome> (between the arms, with
   neighbor(k)), _self_<outcome> (within an arm, with ai(k)). */
mata:
// rows of donors at most as far from score p as the k-th nearest of them
real colvector _s2psm_within(real colvector donors, real colvector ps, real scalar p, real scalar k)
{
    real colvector d, s

    d = abs(ps[donors] :- p)
    s = sort(d, 1)
    return(donors[selectindex(d :<= s[min((k, rows(s)))])])
}

void _s2psm_rule_between(real scalar k, real scalar ate, string scalar outcomes)
{
    real colvector arm, ps, support, used, nn, weight, targets, donors, pick, isdonor
    real matrix Y, M
    string rowvector names
    real scalar a, i, t

    arm = st_data(., "_treated")
    ps = st_data(., "_pscore")
    support = st_data(., "_support")
    names = tokens(outcomes)
    nn = st_data(., "S2REF_nn")
    weight = st_data(., "S2REF_weight")
    used = J(st_nobs(), 1, 0)
    if (cols(names)) {
        Y = st_data(., names)
        M = st_data(., "S2REF_" :+ names)
    }
    for (a = (ate ? 0 : 1); a <= 1; a++) {
        targets = selectindex(support :== 1 :& arm :== a)
        donors = selectindex(support :== 1 :& arm :== 1 - a)
        if (!rows(donors)) continue
        for (i = 1; i <= rows(targets); i++) {
            t = targets[i]
            pick = _s2psm_within(donors, ps, ps[t], k)
            nn[t] = rows(pick)
            used[pick] = used[pick] :+ 1 / rows(pick)
            if (cols(names)) M[t, .] = mean(Y[pick, .])
        }
    }
    // _weight of a donor is its use as a match (missing when unused); the
    // treated keep 1 without ate
    isdonor = selectindex(support :== 1 :& (arm :== 0 :| (ate :& arm :== 1)))
    weight[isdonor] = editvalue(used[isdonor], 0, .)
    st_store(., "S2REF_nn", nn)
    st_store(., "S2REF_weight", weight)
    if (cols(names)) st_store(., "S2REF_" :+ names, M)
}

void _s2psm_rule_self(real scalar k, real scalar altvar, string scalar outcomes)
{
    real colvector arm, ps, support, block, pick
    real matrix Y, S
    real rowvector m, center
    string rowvector names
    real scalar a, i, t, n

    arm = st_data(., "_treated")
    ps = st_data(., "_pscore")
    support = st_data(., "_support")
    names = tokens(outcomes)
    Y = st_data(., names)
    S = st_data(., "S2REF_self_" :+ names)
    for (a = 0; a <= 1; a++) {
        block = selectindex(support :== 1 :& arm :== a)
        if (rows(block) < 2) continue
        for (i = 1; i <= rows(block); i++) {
            t = block[i]
            pick = _s2psm_within(select(block, block :!= t), ps, ps[t], k)
            n = rows(pick)
            m = mean(Y[pick, .])
            if (altvar) {
                // the conditional variance of Abadie et al. (2004, p. 303)
                center = (n :* m + Y[t, .]) :/ (n + 1)
                m = ((Y[t, .] - center):^2 + colsum((Y[pick, .] :- center):^2)) :/ n
            }
            S[t, .] = m
        }
    }
    st_store(., "S2REF_self_" :+ names, S)
}
end

/* Called by _s2psm_cmp with the command line, after psmatch2 ran and its
   results are in the reference. Does nothing unless the call has ties with
   neighbor(k) or ai(k), k > 1, with replacement, on the score. Then:
   - _n1.._nk: the k nearest donors, as psmatch2 lists them without ties;
   - the variables the rule gives (above);
   - the effects and standard errors those variables imply, by psmatch2's
     formulas (its _mktab, the temporary variables in its storage types),
     and their columns of r(table).
   Where the score has no ties the rule gives psmatch2's own matches, so the
   call is still compared with psmatch2 in full. */
capture program drop _s2psm_ties_rule
program define _s2psm_ties_rule
    capture syntax varlist(fv) [if] [in], ties [OUTcome(varlist) Neighbor(integer 1) ai(integer 0) ///
        PScore(varname) ate ALTVariance WARNings NOREPLacement MAHALanobis(string) *]
    if _rc exit
    if "`noreplacement'`mahalanobis'" != "" | (`neighbor' < 2 & `ai' < 2) exit
    if `"`options'"' != "" | (`ai' > 0 & "`ate'" != "") {
        di as error "_s2psm_ties_rule is not written for: `0'"
        exit 198
    }
    if `neighbor' > 1 {
        local noties : subinstr local 0 "ties" "", word
        capture drop _*
        quietly psmatch2 `noties'
        forvalues m = 1/`neighbor' {
            quietly replace S2REF_n`m' = _n`m'
        }
        mata: _s2psm_rule_between(`neighbor', "`ate'" != "", "`outcome'")
    }
    if `ai' > 1 & "`outcome'" != "" {
        mata: _s2psm_rule_self(`ai', "`altvariance'" != "", "`outcome'")
    }
    local col = 0
    tempname att atu ateval seatt seatu seate m0 n0 n1 var0 var1
    foreach v of local outcome {
        quietly summarize `v' if _treated == 1 & _support == 1
        scalar `n1' = r(N)
        scalar `var1' = r(Var)
        scalar `att' = r(mean)
        quietly summarize S2REF_`v' if _treated == 1 & _support == 1, meanonly
        scalar `att' = `att' - r(mean)
        quietly summarize `v' if _treated == 0 & _support == 1
        scalar `n0' = r(N)
        scalar `var0' = r(Var)
        scalar `m0' = r(mean)
        tempvar w shat q
        if `ai' > 0 {
            * Abadie and Imbens (2006), Theorem 7
            quietly gen `w' = max(S2REF_weight, 0)
            quietly gen `shat' = cond("`altvariance'" == "", ///
                (`ai'/(`ai' + 1))*(`v' - S2REF_self_`v')^2, S2REF_self_`v')
            quietly gen `q' = max(0, _treated*(`v' - S2REF_`v' - `att')^2) ///
                + (1 - _treated)*(`w'^2 - `w'/`neighbor')*`shat'
            quietly summarize `q' if _support == 1, meanonly
            scalar `seatt' = sqrt(r(sum))/`n1'
        }
        else {
            quietly gen double `w' = S2REF_weight^2 if _support == 1
            quietly summarize `w' if _treated == 0, meanonly
            scalar `seatt' = sqrt(`var1'/`n1' + `var0'*r(sum)/`n1'^2)
            if "`ate'" != "" {
                quietly summarize S2REF_`v' if _treated == 0 & _support == 1, meanonly
                scalar `atu' = r(mean) - `m0'
                scalar `ateval' = `att'*`n1'/(`n0' + `n1') + `atu'*`n0'/(`n0' + `n1')
                quietly summarize `w' if _treated == 1, meanonly
                scalar `seatu' = sqrt(`var0'/`n0' + `var1'*r(sum)/`n0'^2)
                quietly gen double `q' = (1 + S2REF_weight)^2 if _support == 1
                quietly summarize `q' if _treated == 1, meanonly
                scalar `seate' = `var1'*r(sum)
                quietly summarize `q' if _treated == 0, meanonly
                scalar `seate' = sqrt((`seate' + `var0'*r(sum))/(`n0' + `n1')^2)
            }
        }
        local effects att
        if "`ate'" != "" local effects att atu ate
        foreach e of local effects {
            local b ``e''
            if "`e'" == "ate" local b `ateval'
            local se `se`e''
            scalar S2R_`e' = `b'
            scalar S2R_`e'_`v' = `b'
            scalar S2R_se`e' = `se'
            scalar S2R_se`e'_`v' = `se'
            local ++col
            matrix S2T_ref[1, `col'] = `b'
            matrix S2T_ref[2, `col'] = `se'
            forvalues r = 3/6 {
                matrix S2T_ref[`r', `col'] = .
            }
            if `se' < . & `se' > 0 {
                matrix S2T_ref[3, `col'] = `b'/`se'
                matrix S2T_ref[4, `col'] = 2*normal(-abs(`b'/`se'))
                matrix S2T_ref[5, `col'] = `b' - S2T_ref[8, `col']*`se'
                matrix S2T_ref[6, `col'] = `b' + S2T_ref[8, `col']*`se'
            }
        }
    }
end

capture program drop _s2psm_cmp
program define _s2psm_cmp
    gettoken name 0 : 0, parse(":")
    gettoken colon 0 : 0, parse(":")
    local name = strtrim(`"`name'"')
    local tol 1e-9
    capture drop _*
    capture drop S2REF_*
    capture quietly psmatch2 `0'
    local rc_ref = _rc
    local scalars : r(scalars)
    foreach s of local scalars {
        scalar S2R_`s' = r(`s')
    }
    local macros depvar exog ai_vartype ai_pscorr_method
    foreach m of local macros {
        local M_`m' `"`r(`m')'"'
    }
    tempname T_ct
    capture matrix S2T_ref = r(table)
    local has_table = !_rc
    local vars
    if !`rc_ref' {
        quietly ds _*
        local vars `r(varlist)'
        foreach v of local vars {
            quietly gen double S2REF`v' = `v'
            local lab_`v' : variable label `v'
            local typ_`v' : type `v'
        }
        local char_ref : char _weight[Type]
        _s2psm_ties_rule `0'
    }
    capture drop `vars'
    capture quietly cpsmatch `0'
    local rc_ct = _rc
    local bad
    if `rc_ref' != `rc_ct' local bad rc:`rc_ref'/`rc_ct'
    else if !`rc_ref' {
        foreach s of local scalars {
            local a = scalar(S2R_`s')
            if !(`a' == r(`s') | reldif(`a', r(`s')) <= `tol') local bad `bad' r(`s'):`a'/`=r(`s')'
        }
        foreach m of local macros {
            if `"`r(`m')'"' != `"`M_`m''"' local bad `bad' r(`m')
        }
        if `has_table' {
            capture matrix `T_ct' = r(table)
            if _rc local bad `bad' r(table):absent
            else if "`: rowfullnames `T_ct''" != "`: rowfullnames S2T_ref'" | ///
                "`: colfullnames `T_ct''" != "`: colfullnames S2T_ref'" local bad `bad' r(table):names
            else {
                mata: st_local("differ", strofreal(sum((st_matrix("`T_ct'") :!= st_matrix("S2T_ref")) :& ///
                    !(reldif(st_matrix("`T_ct'"), st_matrix("S2T_ref")) :<= `tol'))))
                if `differ' local bad `bad' r(table):`differ'
            }
        }
        foreach v of local vars {
            capture confirm variable `v'
            if _rc {
                local bad `bad' `v':absent
                continue
            }
            quietly count if !(`v' == S2REF`v' | reldif(`v', S2REF`v') <= `tol')
            if r(N) local bad `bad' `v':`r(N)'
            local lab : variable label `v'
            if subinstr(`"`lab_`v''"', "psmatch2:", "cpsmatch:", 1) != `"`lab'"' local bad `bad' `v':label
            if "`v'" == "_id" | regexm("`v'", "^_n[0-9]+$") {
                local t : type `v'
                if "`t'" != "`typ_`v''" local bad `bad' `v':type
            }
        }
        local char_ct : char _weight[Type]
        if `"`char_ct'"' != `"`char_ref'"' local bad `bad' char
    }
    if "`bad'" == "" & `rc_ref' == 0 test_pass "`name'"
    else if "`bad'" == "" test_pass "`name': same error r(`rc_ref')"
    else test_fail "`name'" "`bad'"
end

* Continuous covariates, a binary covariate and a tied (rounded) score; the
* data in a fixed random order. y2 misses a value (the outcomes' common
* sample), half selects a subsample.
clear
quietly set obs 400
gen double x1 = rnormal()
gen double x2 = rnormal()
gen byte xb = runiform() < .4
gen byte t = (x1 + .5*xb + rnormal() > 0)
gen double y = x1 + t + rnormal()
gen double y2 = x2 + 2*t + rnormal()
quietly replace y2 = . in 7
gen double ps = round(normal(x1/2), .02)
gen byte half = mod(_n, 3) > 0
gen double u = runiform()
sort u

* PSM-12/13/14: nearest neighbor k = 1, 2, 3, with and without ties, on a
* continuous and a tied score: _id, _n#, _nn, _pdif, _<outcome>, labels.
* neighbor(2) ties and neighbor(3) ties are compared with cpsmatch's
* documented rule (_s2psm_ties_rule), which on the tied score is not what
* psmatch2 does.
foreach o in "" "neighbor(2)" "neighbor(3)" "ties" "neighbor(2) ties" "neighbor(3) ties" ///
    "caliper(.005)" "neighbor(3) caliper(.01)" {
    _s2psm_cmp PSM-12 nearest `o', estimated: t x1 x2 xb, outcome(y y2) `o'
    _s2psm_cmp PSM-12 nearest `o', tied score: t, pscore(ps) outcome(y y2) `o'
    _s2psm_cmp PSM-12 nearest `o', if: t x1 x2 if half, outcome(y) `o'
}
_s2psm_cmp PSM-12 nearest, no outcome: t x1 x2
_s2psm_cmp PSM-12 nearest, in range: t x1 x2 in 20/350, outcome(y2) neighbor(2)
* _id outside the sample when the score is given: those rows have a score
* and psmatch2 orders them by it (by its negative with descending)
_s2psm_cmp PSM-12 nearest, given score, if: t if half, pscore(ps) outcome(y)
_s2psm_cmp PSM-12 nearest, given score, in range: t in 20/350, pscore(ps) outcome(y y2) neighbor(2) ate
_s2psm_cmp PSM-12 nearest, given score, if, descending: t if half, pscore(ps) outcome(y) noreplacement descending
* The departure itself, so that a psmatch2 that changes it is noticed (then
* docs/COMPATIBILITY.md and _s2psm_ties_rule need another look): on the tied
* score psmatch2's neighbor(3) ties lists a donor twice among _n1.._n3 for
* some treated rows, and for others stops at 3 donors although more are tied
* at the third distance; cpsmatch lists distinct donors and counts them all.
capture noisily {
    capture drop _*
    quietly psmatch2 t, pscore(ps) outcome(y) neighbor(3) ties
    quietly count if _treated == 1 & (_n1 == _n2 | _n2 == _n3 | _n1 == _n3)
    assert r(N) > 0
    quietly gen long psm_nn = _nn
    capture drop _*
    quietly cpsmatch t, pscore(ps) outcome(y) neighbor(3) ties
    assert _n1 != _n2 & _n2 != _n3 & _n1 != _n3 if _treated == 1
    quietly count if _treated == 1 & psm_nn == 3 & _nn > 3
    assert r(N) > 0
    * every donor of the third-nearest distance or closer, once: the weights
    * of a treated row's donors sum to one, so the controls' sum to the
    * number of treated rows
    quietly summarize _weight if _treated == 0, meanonly
    local used = r(sum)
    quietly count if _treated == 1 & _support == 1
    assert reldif(`used', r(N)) < 1e-12
}
if _rc test_fail "PSM-12 ties with neighbor(3): distinct donors through the third distance" "r(`=_rc')"
else test_pass "PSM-12 ties with neighbor(3): distinct donors through the third distance"
capture drop psm_nn

* PSM-11: without replacement (1-to-1 only, as psmatch2), descending
foreach o in "noreplacement" "noreplacement ties" "noreplacement descending" ///
    "noreplacement descending ties" "noreplacement caliper(.01)" {
    _s2psm_cmp PSM-11 `o', estimated: t x1 x2 xb, outcome(y y2) `o'
    _s2psm_cmp PSM-11 `o', tied score: t, pscore(ps) outcome(y) `o'
}
_s2psm_cmp PSM-11 noreplacement neighbor(2) rejected: t x1 x2, outcome(y) noreplacement neighbor(2)
_s2psm_cmp PSM-11 noreplacement kernel rejected: t x1 x2, outcome(y) noreplacement kernel
_s2psm_cmp PSM-11 noreplacement mahalanobis rejected: t, mahalanobis(x1 x2) outcome(y) noreplacement
_s2psm_cmp PSM-11 descending without noreplacement: t x1 x2, outcome(y) descending
* psmatch2 accepts noreplacement with radius, which ignores it
_s2psm_cmp PSM-11 noreplacement radius: t x1 x2, outcome(y) noreplacement radius caliper(.02)
* Treated outnumbering the controls: without replacement the controls run out
gen byte t2 = (x1 + rnormal() > -1)
_s2psm_cmp PSM-11 controls run out: t2, pscore(ps) outcome(y) noreplacement
_s2psm_cmp PSM-11 controls run out, descending: t2 x1 x2, outcome(y) noreplacement descending

* PSM-02: ate (ATU, ATE, their r() and _<outcome>, _weight, _support)
foreach o in "ate" "ate neighbor(3)" "ate ties neighbor(2)" "ate caliper(.005)" "ate noreplacement" ///
    "ate noreplacement descending" "ate common" "ate kernel" "ate kernel kerneltype(normal)" ///
    "ate radius caliper(.02)" "ate llr" "ate llr kerneltype(biweight)" {
    _s2psm_cmp PSM-02 `o', estimated: t x1 x2 xb, outcome(y y2) `o'
    _s2psm_cmp PSM-02 `o', tied score: t, pscore(ps) outcome(y) `o'
}
_s2psm_cmp PSM-02 ate, many neighbors: t2 x1 x2, outcome(y) ate neighbor(10)
_s2psm_cmp PSM-02 ate, no outcome: t x1 x2, ate

* PSM-16: radius without caliper() (or with caliper(0)) matches every control,
* as psmatch2's missing caliper does
foreach o in "radius" "radius caliper(0)" "radius common" "radius ate" "radius mahalanobis(x1 x2)" {
    _s2psm_cmp PSM-16 `o', estimated: t x1 x2, outcome(y y2) `o'
    _s2psm_cmp PSM-16 `o', tied score: t, pscore(ps) outcome(y) `o'
}

* Radius and kernel matching, each kernel type, and common support
foreach o in "radius caliper(.01)" "radius caliper(.05)" "kernel" "kernel kerneltype(normal)" ///
    "kernel kerneltype(biweight)" "kernel kerneltype(uniform)" "kernel kerneltype(tricube)" ///
    "kernel kerneltype(epan) bwidth(.1)" "kernel kerneltype(normal) bwidth(.02)" ///
    "common" "common neighbor(2)" "common kernel" "common radius caliper(.02)" {
    _s2psm_cmp radius/kernel `o', estimated: t x1 x2 xb, outcome(y y2) `o'
    _s2psm_cmp radius/kernel `o', tied score: t, pscore(ps) outcome(y) `o'
}

* PSM-04: llr, the Epanechnikov kernel (lpoly-smoothed outcomes matched to
* the nearest neighbor, _s_<outcome>) and the other kernels
foreach o in "llr" "llr bwidth(.1)" "llr kerneltype(normal)" "llr kerneltype(biweight)" ///
    "llr kerneltype(uniform) bwidth(.1)" "llr kerneltype(tricube)" "llr neighbor(2)" {
    _s2psm_cmp PSM-04 `o', estimated: t x1 x2 xb, outcome(y y2) `o'
    _s2psm_cmp PSM-04 `o', tied score: t, pscore(ps) outcome(y) `o'
}
_s2psm_cmp PSM-04 llr without outcome: t x1 x2, llr
_s2psm_cmp PSM-04 llr with kernel rejected: t x1 x2, outcome(y) llr kernel

* PSM-03: mahalanobis() (nearest neighbor, kernel, llr, radius) and w()
foreach o in "" "neighbor(3)" "caliper(.5)" "ate" "ate neighbor(2) caliper(.5)" "kernel" ///
    "kernel kerneltype(normal) bwidth(.5)" "kernel ate bwidth(.5)" "llr bwidth(.5)" ///
    "llr ate kerneltype(normal) bwidth(1)" "radius caliper(.3)" "ties" {
    _s2psm_cmp PSM-03 mahalanobis `o': t, mahalanobis(x1 x2) outcome(y y2) `o'
    _s2psm_cmp PSM-03 mahalanobis `o', with score: t x1 x2, mahalanobis(x1 xb) outcome(y) `o'
}
_s2psm_cmp PSM-03 mahalanobis common: t x1, mahalanobis(x1 x2) outcome(y) common
_s2psm_cmp PSM-03 mahalanobis common without score: t, mahalanobis(x1 x2) outcome(y) common
matrix S2W = (2, .5 \ .5, 1)
matrix colnames S2W = x1 x2
matrix rownames S2W = x1 x2
_s2psm_cmp PSM-03 w(): t, mahalanobis(x1 x2) w(S2W) outcome(y) neighbor(2)
_s2psm_cmp PSM-03 w() kernel: t, mahalanobis(x1 x2) w(S2W) outcome(y) kernel bwidth(.5)
_s2psm_cmp PSM-03 w() ate: t, mahalanobis(x1 x2) w(S2W) outcome(y) ate

* PSM-06: ai(#) with samplevar and altvariance (Abadie-Imbens variance),
* on the score and on the Mahalanobis metric. psmatch2's altvariance fails
* with several outcomes, so it has one. On the estimated score the plain
* ai(1), ai(2) and ai(2) ate cases get the Abadie-Imbens (2016) correction.
* With ties, the same-arm matches of ai(2) follow cpsmatch's documented rule
* (_s2psm_ties_rule); "ai(2) altvariance ties" has them on the tied score.
foreach o in "ai(1)" "ai(2)" "ai(2) ate" "ai(2) samplevar" "ai(3) samplevar ate neighbor(2)" ///
    "ai(2) ties" "ai(2) caliper(.01) ate" "ai(2) noreplacement" {
    _s2psm_cmp PSM-06 `o', score: t x1 x2 xb, outcome(y y2) `o'
    _s2psm_cmp PSM-06 `o', mahalanobis: t, mahalanobis(x1 x2) outcome(y y2) `o'
}
foreach o in "ai(2) altvariance" "ai(3) altvariance ate" "ai(2) altvariance samplevar ate" ///
    "ai(2) altvariance ties" {
    _s2psm_cmp PSM-06 `o', score: t, pscore(ps) outcome(y) `o'
    _s2psm_cmp PSM-06 `o', mahalanobis: t, mahalanobis(x1 x2) outcome(y) `o'
}
_s2psm_cmp PSM-06 ai(2) ties, tied score: t, pscore(ps) outcome(y y2) ai(2) ties
* ai(#) above 1 with noreplacement: psmatch2 takes a row's first same-arm
* match without replacement and the others around it, used or not. The
* population variance does not depend on them when every control weighs 1;
* samplevar, ate and altvariance do. With ties psmatch2 can count a
* same-arm match twice, and cpsmatch follows it (no replacement, no rule).
foreach o in "ai(2) noreplacement samplevar" "ai(3) noreplacement ate" "ai(3) noreplacement ate samplevar" ///
    "ai(2) noreplacement descending samplevar" "ai(2) noreplacement caliper(.01) ate samplevar" ///
    "ai(2) noreplacement common samplevar" "ai(2) noreplacement ties" "ai(3) noreplacement ties ate samplevar" {
    _s2psm_cmp PSM-06 `o', score: t x1 x2 xb, outcome(y y2) `o'
    _s2psm_cmp PSM-06 `o', tied score: t, pscore(ps) outcome(y y2) `o'
}
_s2psm_cmp PSM-06 ai(1) noreplacement samplevar: t x1 x2 xb, outcome(y y2) ai(1) noreplacement samplevar
_s2psm_cmp PSM-06 ai(2) noreplacement altvariance: t x1 x2 xb, outcome(y) ai(2) noreplacement altvariance
_s2psm_cmp PSM-06 ai(3) noreplacement altvariance ties, tied score: t, pscore(ps) outcome(y) ai(3) noreplacement altvariance ties ate
_s2psm_cmp PSM-06 ai(2) noreplacement, controls run out: t2, pscore(ps) outcome(y) ai(2) noreplacement samplevar
* That pass is the plugin's: it says so (_cpsmatch_ai_pass, without which
* cpsmatch stops with r(601)), and every row gets its same-arm matches
capture noisily {
    capture drop _*
    quietly cpsmatch t x1 x2 xb, outcome(y) ai(3) noreplacement samplevar
    assert scalar(_cpsmatch_ai_pass) == 2
    assert _self_y != 0 & _self_y < . if _support == 1
}
if _rc test_fail "PSM-06 ai(3) noreplacement: the plugin's same-arm pass" "r(`=_rc')"
else test_pass "PSM-06 ai(3) noreplacement: the plugin's same-arm pass"

* ai() when an arm has one row on support. psmatch2's same-arm pass then
* reads the row after that one in its order and takes it as the row's match
* if it has a score within the caliper: the first row of the other arm
* (one1c: one control), a row with another treatment code (one1t), a row off
* support (common ate) or one outside the sample (if, with a given score;
* its y2 can be missing). cpsmatch does the same.
gen byte one1c = _n != 1
gen byte one1t = cond(_n == 1, 1, cond(mod(_n, 7) == 0, 2, 0))
gen byte last1t = _n == 1
foreach o in "ai(1)" "ai(2) ate samplevar" "ai(2) altvariance ate" "ai(2) ties" ///
    "ai(2) noreplacement samplevar" "ai(3) noreplacement descending ate samplevar" ///
    "ai(2) caliper(.3) ate samplevar" "ai(2) common samplevar" {
    _s2psm_cmp PSM-06 one control, `o': one1c, pscore(ps) outcome(y) `o'
    _s2psm_cmp PSM-06 one treated before another code, `o': one1t, pscore(ps) outcome(y) `o'
}
_s2psm_cmp PSM-06 one control, correction: one1c x1 x2, outcome(y y2) ai(2) ate debug
_s2psm_cmp PSM-06 one control, no outcome: one1c, pscore(ps) ai(2)
_s2psm_cmp PSM-06 one treated before a row off support: last1t, pscore(ps) outcome(y) ai(2) samplevar common ate
_s2psm_cmp PSM-06 one treated before a row outside the sample: last1t if half | _n == 1, pscore(ps) outcome(y) ai(2) samplevar
_s2psm_cmp PSM-06 one treated before a row without y2: last1t, pscore(ps) outcome(y y2) ai(2) samplevar
_s2psm_cmp PSM-06 one treated before a row without a score: last1t x1 x2 if half | _n == 1, outcome(y) ai(2) samplevar
* When no row follows (the treated row is the data's last in psmatch2's
* order), and with mahalanobis() and w(), psmatch2 stops with a Mata error,
* r(3301); cpsmatch stops with r(2001) and a message (docs/COMPATIBILITY.md).
* Without w() both stop with r(2001).
capture noisily {
    foreach c in "last1t, pscore(ps) outcome(y) ai(2)" "last1t x1 x2, outcome(y) ai(2) samplevar" ///
        "last1t, pscore(ps) outcome(y) ai(2) noreplacement ate" "last1t, pscore(ps) ai(2)" ///
        "last1t, mahalanobis(x1 x2) w(S2W) outcome(y) ai(2)" ///
        "one1c, mahalanobis(x1 x2) w(S2W) outcome(y) ai(2) ate" {
        capture drop _*
        capture psmatch2 `c'
        assert _rc == 3301
        capture drop _*
        capture cpsmatch `c'
        assert _rc == 2001
    }
}
if _rc test_fail "PSM-06 ai with a one-row arm and no row to borrow: r(2001)" "r(`=_rc')"
else test_pass "PSM-06 ai with a one-row arm and no row to borrow: r(2001)"
_s2psm_cmp PSM-06 one control, mahalanobis: one1c, mahalanobis(x1 x2) outcome(y) ai(2) ate
_s2psm_cmp PSM-06 one treated, mahalanobis: last1t, mahalanobis(x1 x2) outcome(y) ai(2)
_s2psm_cmp PSM-06 ai with kernel rejected: t x1 x2, outcome(y) ai(2) kernel
_s2psm_cmp PSM-06 ai with radius: t x1 x2, outcome(y) ai(2) radius caliper(.1)
* psmatch2 4.1.0 has samplevar where 4.0.12 had population
_s2psm_cmp PSM-06 population rejected: t x1 x2, outcome(y) ai(2) population
_s2psm_cmp PSM-06 samplevar without ai: t x1 x2, outcome(y) samplevar
_s2psm_cmp PSM-06 samplevar, given score: t, pscore(ps) outcome(y) ai(2) samplevar ate

* PSM-17: the Abadie-Imbens (2016) correction for an estimated score
* (r(ai_pscorr) is 1): probit and logit, several neighbors, with and without
* ate (without it the controls' _nn is missing and their _pdif filled),
* debug's parts of the correction, factor variables and an omitted
* covariate in the score model, a subsample, several outcomes
foreach o in "ai(2) debug" "ai(2) ate debug" "ai(1) logit" "ai(3) neighbor(3)" ///
    "ai(3) neighbor(3) ate logit debug" "ai(2) neighbor(2) quietly" "ai(2) warnings" {
    _s2psm_cmp PSM-17 correction `o': t x1 x2 xb, outcome(y y2) `o'
}
_s2psm_cmp PSM-17 correction, factor variables: t x1 c.x2##i.xb, outcome(y) ai(2) ate debug
_s2psm_cmp PSM-17 correction, if: t x1 x2 if half, outcome(y) ai(2) debug
_s2psm_cmp PSM-17 correction, in range: t x1 x2 xb in 20/350, outcome(y2) ai(2) ate
_s2psm_cmp PSM-17 correction, no outcome: t x1 x2, ai(2)
_s2psm_cmp PSM-17 correction, no outcome, ate: t x1 x2, ai(2) ate
gen double x1twice = 2*x1
_s2psm_cmp PSM-17 correction, omitted covariate: t x1 x1twice x2, outcome(y) ai(2) ate debug
drop x1twice
* The options that leave the correction out (r(ai_pscorr) is 0), and debug
* and r(ai) without ai()
foreach o in "ai(2) index" "ai(2) odds logit" "ai(2) common" "ai(2) trim(5)" "ai(2) samplevar debug" ///
    "ai(2) caliper(.05) debug" "debug" "debug ate kernel" {
    _s2psm_cmp PSM-17 no correction `o': t x1 x2 xb, outcome(y y2) `o'
}
_s2psm_cmp PSM-17 no correction, given score: t, pscore(ps) outcome(y) ai(2) ate debug
* r(table) follows the confidence level
set level 90
_s2psm_cmp PSM-17 r(table) at level 90: t x1 x2 xb, outcome(y y2) ate
_s2psm_cmp PSM-17 r(table) at level 90, correction: t x1 x2 xb, outcome(y) ai(2) ate
set level 95

* PSM-07: trim(#), with ate and with the other methods
foreach o in "trim(5)" "trim(0)" "trim(20) ate" "trim(10) kernel" "trim(10) common neighbor(2)" ///
    "trim(5) ate ai(2)" {
    _s2psm_cmp PSM-07 `o', estimated: t x1 x2 xb, outcome(y y2) `o'
    _s2psm_cmp PSM-07 `o', given score: t, pscore(ps) outcome(y) `o'
}
_s2psm_cmp PSM-07 trim out of range: t, pscore(ps) outcome(y) trim(-5)
_s2psm_cmp PSM-07 trim above 100 with common: t, pscore(ps) outcome(y) trim(150) common
_s2psm_cmp PSM-07 trim without a score: t, mahalanobis(x1) outcome(y) trim(5)

* S3: odds, index, quietly, warnings, logit
foreach o in "odds" "odds logit" "odds logit trim(5) ate" "index" "index common" "quietly" ///
    "quietly logit kernel" "warnings" "warnings ate" "warnings neighbor(2) ties" {
    _s2psm_cmp S3 `o': t x1 x2 xb, outcome(y y2) `o'
}
_s2psm_cmp S3 warnings, tied score: t, pscore(ps) outcome(y) warnings
_s2psm_cmp S3 warnings, tied score, ate: t, pscore(ps) outcome(y) warnings ate

* PSM-05: spline needs the -spline- package, as in psmatch2
capture which spline
if _rc {
    _s2psm_cmp PSM-05 spline without the spline package: t x1 x2, outcome(y) spline
}
else {
    _s2psm_cmp PSM-05 spline: t x1 x2, outcome(y) spline
    _s2psm_cmp PSM-05 spline nknots(): t x1 x2, outcome(y) spline nknots(3) ate
}

* Other treatment codes and one-arm samples
gen byte t3 = cond(mod(_n, 5) == 0, 2, t)
foreach o in "" "ate" "kernel ate" "ai(2)" "neighbor(2) ties" {
    _s2psm_cmp other codes `o': t3, pscore(ps) outcome(y) `o'
}
_s2psm_cmp other codes, mahalanobis ate: t3, mahalanobis(x1 x2) outcome(y) ate
gen byte t0 = 0
gen byte t1 = 1
_s2psm_cmp no treated, ate: t0, pscore(ps) outcome(y) ate
_s2psm_cmp no controls, ate kernel: t1, pscore(ps) outcome(y) ate kernel

* Option checks in psmatch2's order
_s2psm_cmp checks: no score: t, outcome(y)
_s2psm_cmp checks: score and covariates: t x1, pscore(ps) outcome(y)
_s2psm_cmp checks: kerneltype: t, pscore(ps) outcome(y) kernel kerneltype(foo)
_s2psm_cmp checks: neighbor(0): t, pscore(ps) outcome(y) neighbor(0)

* Contracts: value labels, _weight's characteristic, _nn on support
capture noisily {
    capture drop _*
    quietly cpsmatch t x1 x2, outcome(y) ate
    assert "`: value label _treated'" == "_treated" & "`: value label _support'" == "_support"
    assert "`: label _treated 0'" == "Untreated" & "`: label _treated 1'" == "Treated"
    assert "`: label _support 0'" == "Off support" & "`: label _support 1'" == "On support"
    assert "`: char _weight[Type]'" == "aweight"
    assert _nn < . if _support == 1
    assert !missing(r(atu), r(ate), r(atu_y), r(ate_y))
    assert !missing(r(seatu), r(seate), r(seatu_y), r(seate_y))
    assert r(ai) == 0 & "`r(ai_vartype)'" == ""
    assert rowsof(r(table)) == 9 & colsof(r(table)) == 3
    quietly cpsmatch t x1 x2, outcome(y) ai(2) ate
    assert !missing(r(seatu), r(seate), r(seatu_y), r(seate_y))
    assert r(ai) == 2 & r(ai_pscorr) == 1 & "`r(ai_vartype)'" == "population"
    assert r(method) == "nearest" & r(metric) == "pscore"
    quietly cpsmatch t x1 x2, outcome(y) ai(2) ate samplevar
    assert r(ai_pscorr) == 0 & "`r(ai_vartype)'" == "sample"
}
if _rc test_fail "PSM-14 labels, characteristic and r()" "r(`=_rc')"
else test_pass "PSM-14 labels, characteristic and r()"

* Several outcomes: each outcome's matched values and effects equal a
* run with that outcome alone, over the common sample
capture noisily {
    capture drop _*
    quietly cpsmatch t x1 x2 if !missing(y2), outcome(y) ai(2) ate
    local a = r(ate)
    local s = r(seate)
    quietly gen double one_y = _y
    capture drop _*
    quietly cpsmatch t x1 x2, outcome(y y2) ai(2) ate
    assert reldif(r(ate_y), `a') < 1e-12 & reldif(r(seate_y), `s') < 1e-12
    assert _y == one_y | (missing(_y) & missing(one_y))
    drop one_y
}
if _rc test_fail "PSM-01 several outcomes with ate and ai()" "r(`=_rc')"
else test_pass "PSM-01 several outcomes with ate and ai()"

* Several outcomes when matching takes rows off the support that common or
* trim() left (a caliper, controls running out): the match is the same for
* every outcome, from the support before matching (_id, _n#, _nn)
foreach o in "common caliper(.005)" "common caliper(.005) ate" "trim(10) caliper(.005)" ///
    "trim(10) caliper(.005) ate ai(2) samplevar" "common caliper(.005) noreplacement" ///
    "common caliper(.01) ai(2) noreplacement samplevar" {
    _s2psm_cmp PSM-01 several outcomes, `o': t x1 x2 xb, outcome(y y2) `o'
    _s2psm_cmp PSM-01 several outcomes, `o', tied score: t, pscore(ps) outcome(y y2) `o'
}
_s2psm_cmp PSM-01 several outcomes, common, controls run out: t2, pscore(ps) outcome(y y2) common noreplacement

capture drop _*
capture drop S2REF_*
print_summary s2_psm
if $TESTS_FAILED > 0 exit 9
global CTOOLS_COMPONENT_COMPLETE "s2_psm"
    global CTOOLS_CASE_COMPLETE "s2_psm"
    exit
}

di as error "Unknown s2_psm case: `__case'"
exit 198
