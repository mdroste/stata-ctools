*! version 1.2.6 20261111 github.com/mdroste/stata-ctools
*! The wild bootstrap table after the coefficient table of an HDFE estimator.
* Design: docs/development/BOOTSTRAP_PLAN.md (section 4.7).
*
* _ctools_wild_table
*     Shows, from the current e(): the tested coefficients with e(b), their
*     analytic t (e(wb_t)), the bootstrap p-values (e(wb_p)) and, when
*     posted, the confidence sets by test inversion (e(wb_ci),
*     e(wb_cistatus)); for the score bootstrap the score statistic
*     (e(wb_score)) in place of e(wb_t); then notes on the p-value type, enumeration, few
*     Rademacher clusters, open, disjoint or empty sets, trial values whose p could
*     not be evaluated, restricted fits that failed
*     (e(wb_fail)) and replicates without a statistic (e(wb_misreps)).
*     e(wb_method): wcr, wcu (creghdfe), wre, wue (civreghdfe) or score
*     (cppmlhdfe).

program define _ctools_wild_table
    version 14.1
    tempname t p ci st fail
    matrix `t' = e(wb_t)
    * the score bootstrap studentizes the score statistic, not the Wald z
    local score_t = 0
    if "`e(wb_method)'" == "score" {
        capture confirm matrix e(wb_score)
        if !_rc {
            matrix `t' = e(wb_score)
            local score_t = 1
        }
    }
    matrix `p' = e(wb_p)
    local hasci = 0
    capture confirm matrix e(wb_ci)
    if !_rc {
        local hasci = 1
        matrix `ci' = e(wb_ci)
        matrix `st' = e(wb_cistatus)
    }
    local hasfail = 0
    capture confirm matrix e(wb_fail)
    if !_rc {
        local hasfail = 1
        matrix `fail' = e(wb_fail)
    }
    if "`e(wb_method)'" == "score" {
        local kind = cond("`e(wb_bootclust)'" != "", "Score cluster bootstrap", "Score bootstrap")
        local how null imposed
    }
    else {
        local kind = cond("`e(wb_bootclust)'" != "", "Wild cluster bootstrap", "Wild bootstrap")
        local how = cond("`e(wb_method)'" == "wcr", "restricted (WCR)", ///
            cond("`e(wb_method)'" == "wre", "restricted efficient (WRE)", ///
            cond("`e(wb_method)'" == "wue", "unrestricted efficient (WUE)", "unrestricted (WCU)")))
    }
    local wt `e(wb_weighttype)'
    if inlist("`wt'", "rademacher", "mammen", "webb") local wt = strproper("`wt'")
    local who = cond("`e(wb_bootclust)'" != "", "(`e(wb_bootclust)')", "(observations)")
    local stat = cond("`e(wb_method)'" == "score", "z", "t")
    local pl = cond("`e(wb_ptype)'" == "lower", "P<`stat'", cond("`e(wb_ptype)'" == "upper", "P>`stat'", "P>|`stat'|"))
    di as text ""
    di as text "`kind': `how', `wt' weights"
    di as text "Replications = " as result strtrim(string(e(wb_reps), "%12.0fc")) as text ///
        _col(40) "Bootstrap clusters = " as result strtrim(string(e(wb_N_clust), "%12.0fc")) ///
        as text " `who'"
    di as text "{hline 13}{c TT}{hline 64}"
    di as text %12s abbrev("`e(depvar)'", 12) " {c |}" _col(20) "Coef." _col(34) "`stat'" _col(39) "`pl'" ///
        _col(50) "[`=strofreal(e(wb_level))'% Conf. Interval]"
    di as text "{hline 13}{c +}{hline 64}"
    local names : colnames `p'
    local j 0
    local open
    local disjoint
    local empty
    local trialfail
    local failed
    foreach c of local names {
        local ++j
        local lo .
        local hi .
        if `hasci' {
            local lo = `ci'[`j', 1]
            local hi = `ci'[`j', 2]
            if `st'[`j', 1] | `st'[`j', 2] local open `open' `c'
            if `st'[`j', 3] local disjoint `disjoint' `c'
            if `st'[`j', 4] local empty `empty' `c'
            if colsof(`st') >= 5 {
                if `st'[`j', 5] local trialfail `trialfail' `c'
            }
        }
        if `hasfail' {
            if `fail'[1, `j'] local failed `failed' `c'
        }
        di as text %12s abbrev("`c'", 12) " {c |}" as result _col(15) %10.0g _b[`c'] ///
            _col(27) %8.2f `t'[1, `j'] _col(38) %6.4f `p'[1, `j'] _col(49) %10.0g `lo' ///
            _col(62) %10.0g `hi'
    }
    di as text "{hline 13}{c BT}{hline 64}"
    local pdesc = cond("`e(wb_ptype)'" == "symmetric", "symmetric", cond("`e(wb_ptype)'" == "equaltail", "equal-tailed", "`e(wb_ptype)'-tailed"))
    di as text "P-values are `pdesc' bootstrap-`stat'" _c
    if `hasci' di as text "; intervals invert the bootstrap test." _c
    di as text ""
    if `score_t' di as text "z is the score statistic (e(wb_score)); e(wb_t) holds the Wald z."
    if e(wb_enumerated) == 1 {
        di as text "All " as result e(wb_reps) as text " Rademacher sign vectors were enumerated (no random draws)."
    }
    if "`e(wb_weighttype)'" == "rademacher" & e(wb_N_clust) <= 11 {
        di as text "{p 0 6 2}Note: with " e(wb_N_clust) " bootstrap clusters, Rademacher weights allow only " ///
            2^e(wb_N_clust) " distinct draws; consider weighttype(webb).{p_end}"
    }
    if "`e(wb_fe_path)'" == "partialled" {
        di as text "{p 0 6 2}Note: absorbed effects not nested in the bootstrap clusters enter the bootstrap as partialled out (fedgp(partialled)).{p_end}"
    }
    if "`open'" != "" di as text "{p 0 6 2}Note: the confidence set of`open' is unbounded (missing bound) within 1000 standard errors.{p_end}"
    if "`disjoint'" != "" di as text "{p 0 6 2}Note: the confidence set of`disjoint' is not an interval; the bounds enclose it.{p_end}"
    if "`empty'" != "" di as text "{p 0 6 2}Note: no trial value of`empty' was accepted, or its interval could not be searched: no bounds.{p_end}"
    if "`trialfail'" != "" di as text "{p 0 6 2}Note: the p-value could not be evaluated at some trial values for`trialfail'; they were treated as not rejected, so the interval may be wider than the confidence set.{p_end}"
    if "`failed'" != "" di as text "{p 0 6 2}Note: the fit under the null did not converge for`failed'; no bootstrap result.{p_end}"
    if e(wb_misreps) > 0 & e(wb_misreps) < . {
        di as text "Note: " e(wb_misreps) " replicates had no defined statistic and were left out."
    }
end
