*! version 1.3.2 20261124 github.com/mdroste/stata-ctools
*! original: psmatch2 4.1.0 (June 11, 2026) by Edwin Leuven and Barbara Sianesi
*! original: see help cpsmatch for details

program define cpsmatch, rclass sortpreserve
    version 14.1

    * Check observation limit (Stata plugin API limitation)
    if _N > 2147483647 {
        di as error "ctools does not support datasets exceeding 2^31 (2.147 billion) observations"
        di as error "This is a limitation of Stata's plugin API"
        exit 920
    }

    capture noisily syntax varlist(min=1 numeric fv) [if] [in], ///
        [OUTcome(varlist numeric) Pscore(varname numeric) ///
         Logit Probit INDEX ODDS QUIetly ///
         Neighbor(integer 1) Caliper(real 0) Radius Kernel LLR ///
         Kerneltype(string) Bwidth(real 0.06) SPLINE NKnots(integer 0) ///
         MAHALanobis(varlist numeric fv) W(string) ///
         Common TRIM(real 100) NOREPlacement Ties Descending ///
         ATE AI(integer 0) SAMPLEVar ALTVariance DEBUG WARNings ///
         Verbose THReads(integer 0)]
    if _rc == 109 {
        exit 2000
    }
    else if _rc != 0 {
        exit _rc
    }
    * The if is evaluated now, before any helper replaces the r() or e()
    * results it may name; the matching sample is completed below.
    marksample touse

    * =========================================================================
    * UPFRONT VALIDATION
    * =========================================================================

    * Split varlist into treatment variable and covariates
    gettoken depvar treatvar : varlist
    * depvar = treatment variable (first var)
    * treatvar = covariates (remaining vars); factor variables are allowed,
    * as in psmatch2, and expanded by the score model
    _fv_check_depvar `depvar'
    local treatvar = strtrim(`"`treatvar'"')
    local k_covariates : word count `treatvar'

    * psmatch2's checks, in its order. The method is settled in two steps
    * as there: radius counts as nearest-neighbor matching for the
    * noreplacement and ai() checks, and only then becomes kernel matching.
    local metric = cond("`mahalanobis'" != "", "mahalanobis", "pscore")
    if `k_covariates' == 0 & "`pscore'" == "" & "`metric'" == "pscore" {
        di as error "You should either specify a " as input "varlist" as error " or " as input "propensity score"
        exit 198
    }
    if `k_covariates' > 0 & "`pscore'" != "" {
        di as error "You cannot specify both a " as input "varlist" as error " AND a " as input "propensity score"
        exit 198
    }
    if "`ties'" != "" & "`metric'" == "mahalanobis" {
        di as text "Ties is not implemented for Mahalanobis matching. This option will be ignored."
        local ties
    }
    local method "neighbor"
    if "`kernel'" != "" local method "kernel"
    if "`llr'" != "" local method "llr"
    if "`llr'" != "" & "`kernel'" != "" {
        di as error "You cannot do kernel and llr matching at the same time"
        exit 198
    }
    if "`noreplacement'" != "" & ("`method'" != "neighbor" | "`metric'" != "pscore" | `neighbor' > 1) {
        di as error "Matching without replacement is only implemented with 1-to-1 propensity score matching"
        exit 198
    }
    if "`descending'" != "" & "`noreplacement'" == "" {
        di as error "Option " as input "descending" as error " makes only sense when matching without replacement"
        exit 198
    }
    if `ai' > 0 & "`method'" != "neighbor" {
        di as error "Option ai() is only allowed with nearest-neighbor matching"
        exit 198
    }
    * psmatch2 corrects the ai() standard errors for the estimation of the
    * score (Abadie and Imbens 2016) when the score is estimated here and
    * matched on by nearest neighbors, with the population variance and none
    * of the options below. It then matches the controls to the treated as
    * well, as with ate, whether or not ate is specified (work_ate).
    local ai_given = `ai'
    local pscorr = (`ai' > 0 & "`samplevar'" == "" & "`method'" == "neighbor" ///
        & "`metric'" == "pscore" & `k_covariates' > 0 & "`pscore'" == "" ///
        & "`index'`odds'`ties'`noreplacement'`altvariance'`common'" == "" ///
        & `caliper' <= 0 & `trim' == 100 & "`radius'`spline'" == "")
    local work_ate = cond("`ate'" != "" | `pscorr', "ate", "")
    local psfixnote = (`ai' > 0 & "`metric'" == "pscore" & !`pscorr')
    if inlist("`method'", "kernel", "llr") & "`kerneltype'" == "" local kerneltype "epan"
    if !inlist("`kerneltype'", "", "epan", "normal", "biweight", "uniform", "tricube") {
        di as error "Kerneltype `kerneltype' not recognized"
        exit 198
    }
    * psmatch2's radius is kernel matching with the uniform kernel and the
    * caliper as bandwidth. Without caliper() (or with 0) psmatch2's caliper
    * is missing, so every control is within it: a radius no score
    * difference reaches stands for that (and the weights are adjusted after
    * matching, as psmatch2's missing bandwidth spreads them).
    local radius_all = 0
    if "`radius'" != "" {
        local method "kernel"
        local kerneltype "uniform"
        if `caliper' == 0 {
            local caliper = 1e300
            local radius_all = 1
        }
        local bwidth = `caliper'
        * psmatch2 reaches its table without the ai() variables
        if `ai' > 0 & "`outcome'" != "" {
            di as error "_self_`: word 1 of `outcome'' not found"
            exit 111
        }
        local ai = 0
    }
    if "`method'" == "llr" & "`kerneltype'" == "epan" & "`metric'" == "pscore" & "`outcome'" == "" {
        di as error "varlist required"
        exit 100
    }
    * psmatch2's nearest-neighbor matching needs _n1 (r(111) without neighbors)
    if `neighbor' <= 0 & ("`method'" == "neighbor" | ("`method'" == "llr" & "`kerneltype'" == "epan" ///
        & "`metric'" == "pscore") | "`spline'" != "") {
        di as error "variable _n1 not found"
        exit 111
    }

    local kernelcode = 0  /* 0=epan, 1=normal, 2=biweight, 3=uniform, 4=tricube */
    if "`kerneltype'" == "normal" local kernelcode = 1
    else if "`kerneltype'" == "biweight" local kernelcode = 2
    else if "`kerneltype'" == "uniform" local kernelcode = 3
    else if "`kerneltype'" == "tricube" local kernelcode = 4

    * Resolve FV base names and every output before deleting any variables.
    fvrevar `varlist' `mahalanobis', list
    local inputs `r(varlist)' `pscore' `outcome'
    local generated _treated _support _weight _pscore _id _nn _pdif
    forvalues n = 1/`neighbor' {
        local generated `generated' _n`n'
    }
    local willsmooth = "`spline'" != "" | ("`method'" == "llr" & ///
        "`kerneltype'" == "epan" & "`metric'" == "pscore")
    foreach v of local outcome {
        local generated `generated' _`v'
        if `willsmooth' local generated `generated' _s_`v'
        if `ai' > 0 local generated `generated' _self_`v'
    }
    local seen
    foreach v of local generated {
        confirm name `v'
        if `: list v in inputs' {
            di as error "input variable `v' conflicts with a matching output"
            exit 111
        }
        if `: list v in seen' {
            di as error "matching outputs have conflicting name `v'"
            exit 110
        }
        local seen `seen' `v'
    }

    * Verbose timing
    local __do_timing = ("`verbose'" != "")
    if `__do_timing' {
        local __ct_elapsed90 = 0
        local __ct_elapsed91 = 0
        local __ct_elapsed92 = 0
        local __ct_elapsed94 = 0
        quietly _ctools_clock __ct_start90
        quietly _ctools_clock __ct_start91
    }

    * =========================================================================
    * LOAD PLUGIN
    * =========================================================================

    _ctools_load
    * Stata scopes plugin registrations to the calling ado program.
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    capture confirm number 0

    * =========================================================================
    * SAMPLE AND PROPENSITY SCORE
    * =========================================================================

    * The Mahalanobis metric and local linear regression other than the
    * Epanechnikov kernel on the score are matched in Mata, over the
    * observations in their order now; everything else goes to the plugin.
    local __mata = ("`metric'" != "pscore") | ("`method'" == "llr" & "`spline'" == "" ///
        & "`kerneltype'" != "epan")
    tempvar order
    if `__mata' quietly gen long `order' = _n

    * psmatch2's generated variables, dropped in one command: each drop of
    * a variable that is not the last one rewrites the dataset
    local __cand _treated _support _weight _pscore _id _nn _pdif
    forvalues n = 1/`neighbor' {
        local __cand `__cand' _n`n'
    }
    foreach v of local outcome {
        local __cand `__cand' _`v'
    }
    local __old
    foreach v of local __cand {
        capture confirm variable `v', exact
        if !_rc local __old `__old' `v'
    }
    local __old : list uniq __old
    if "`__old'" != "" drop `__old'

    * Matching sample as psmatch2 builds it: treatment, covariates, the
    * outcomes and the Mahalanobis variables nonmissing, restricted to the
    * score model's e(sample), then a nonmissing score.
    markout `touse' `outcome'
    if "`mahalanobis'" != "" {
        fvrevar `mahalanobis', list
        markout `touse' `r(varlist)'
    }

    local score
    if `k_covariates' > 0 {
        local model = cond("`logit'" != "", "logit", "probit")
        quietly count if `touse'
        local __n_model = r(N)
        `quietly' `model' `depvar' `treatvar' if `touse', nolog
        * e(sample) is within the sample: equal when the model used all of it
        if e(N) != `__n_model' quietly replace `touse' = e(sample)
        if "`odds'" != "" | "`common'" != "" | `trim' < 100 {
            * common support uses the score as estimated, before odds
            tempvar score
            quietly predict double `score', `index'
            quietly gen double _pscore = `score' if `touse'
        }
        else {
            quietly predict double _pscore if `touse', `index'
            local score _pscore
        }
        label variable _pscore "cpsmatch: Propensity Score"
        if `pscorr' {
            * The first stage of the correction: the variance of the
            * coefficients, the model's columns (those not omitted; an
            * omitted one has no variance and so no part in it), and the
            * derivative of the score with respect to the index.
            tempname psV
            tempvar psdens
            matrix `psV' = e(V)
            local psnames : colnames e(b)
            if "`model'" == "logit" quietly gen double `psdens' = _pscore*(1 - _pscore) if `touse'
            else {
                tempvar psxb
                quietly predict double `psxb' if `touse', xb
                quietly gen double `psdens' = normalden(`psxb') if `touse'
                drop `psxb'
            }
            local psx
            local pskeep
            local j = 0
            foreach name of local psnames {
                local ++j
                if "`name'" == "_cons" {
                    local pscons = `j'
                    continue
                }
                if `psV'[`j', `j'] == 0 continue
                fvrevar `name'
                local psx `psx' `r(varlist)'
                local pskeep `pskeep' `j'
            }
            * the constant last, as the correction orders its columns
            if "`pscons'" != "" local pskeep `pskeep' `pscons'
            else local pscons = 0
        }
    }
    else if "`pscore'" != "" {
        local score `pscore'
        quietly gen double _pscore = `pscore'
        label variable _pscore "cpsmatch: Propensity Score"
    }
    if "`score'" != "" markout `touse' _pscore
    * odds: match on the log odds; common support still uses the score. A
    * row without log odds (a score of 0 or 1, or with index one outside
    * (0, 1)) leaves the sample, where psmatch2 keeps it unmatchable.
    if "`odds'" != "" & "`score'" != "" {
        quietly replace _pscore = ln(_pscore/(1 - _pscore))
        markout `touse' _pscore
    }

    quietly gen byte _treated = `depvar' if `touse'
    label variable _treated "cpsmatch: Treatment assignment"
    capture label drop _treated
    label define _treated 0 "Untreated" 1 "Treated"
    label values _treated _treated

    * =========================================================================
    * COMMON SUPPORT
    * =========================================================================

    if ("`common'" != "" | `trim' < 100) & "`score'" == "" {
        di as error "With option 'common' a propensity score is needed. Provide one"
        di as error "with option 'pscore()' or estimate one. See the help file for more details."
        exit 198
    }
    quietly gen byte _support = 1 if `touse'
    label variable _support "cpsmatch: Common support"
    capture label drop _support
    label define _support 0 "Off support" 1 "On support"
    label values _support _support
    if "`common'" != "" | `trim' < 100 {
        if !inrange(`trim', 0, 100) {
            di as error "Trim level out of range"
            exit 198
        }
        quietly _cpsmatch_support `score', level(`trim') `ate'
    }

    quietly count if `touse'
    local nobs = r(N)
    if `nobs' == 0 {
        di as error "cpsmatch: no observations"
        exit 2000
    }

    * =========================================================================
    * SMOOTHED OUTCOMES (llr with the Epanechnikov kernel, spline)
    * =========================================================================

    * psmatch2 matches the nearest neighbors' smoothed control outcomes: a
    * local linear (lpoly) or spline fit on the score, over the supported
    * controls (and, with ate, the treated for the treated).
    local smooth = 0
    local matchvars `outcome'
    if "`method'" == "llr" & "`kerneltype'" == "epan" & "`metric'" == "pscore" {
        local method "neighbor"
        local smooth = 1
        local matchvars
        foreach v of local outcome {
            capture drop _s_`v'
            quietly lpoly `v' _pscore if _treated == 0 & _support == 1, nograph degree(1) ///
                at(_pscore) generate(_s_`v') bwidth(`bwidth')
            if "`ate'" != "" {
                tempvar s
                quietly lpoly `v' _pscore if _treated == 1 & _support == 1, nograph degree(1) ///
                    at(_pscore) generate(`s') bwidth(`bwidth')
                quietly replace _s_`v' = `s' if _treated == 1 & _support == 1
                drop `s'
            }
            label variable _s_`v' "cpsmatch: smoothed outcome variable"
            local matchvars `matchvars' _s_`v'
        }
    }
    if "`spline'" != "" {
        capture which spline
        if _rc {
            di as error "You need to install package -spline-, for option spline"
            exit 198
        }
        local method "neighbor"
        local smooth = 2
        if `nknots' == 0 {
            quietly count if _treated == 0 & _support == 1
            local nknots = int(r(N)^0.25)
        }
        local matchvars
        foreach v of local outcome {
            capture drop _s_`v'
            quietly spline `v' _pscore if _treated == 0 & _support == 1, gen(_s_`v') ///
                nknots(`nknots') nograph
            if "`ate'" != "" {
                tempvar s
                quietly spline `v' _pscore if _treated == 1 & _support == 1, gen(`s') ///
                    nknots(`nknots') nograph
                quietly replace _s_`v' = `s' if _treated == 1 & _support == 1
                drop `s'
            }
            label variable _s_`v' "cpsmatch: smoothed outcome variable using -spline-"
            local matchvars `matchvars' _s_`v'
        }
    }

    * =========================================================================
    * CREATE OUTPUT VARIABLES
    * =========================================================================

    * The plugin matches on the propensity score; Mata matches on the
    * Mahalanobis metric and does local linear regression matching.
    local engine = cond("`metric'" == "pscore" & inlist("`method'", "neighbor", "kernel"), "plugin", "mata")
    local nearest = ("`method'" == "neighbor")

    * As psmatch2: the treated weigh 1 (0 with ate, then their use as
    * matches), the matched outcome starts at 0 for the rows to be matched.
    * Mata accumulates into those values. The plugin computes them itself and
    * stores every output on every row of the sample, so its outputs start
    * as missing.
    local __plugin = ("`engine'" == "plugin")
    if "`engine'" == "mata" & !`__mata' {
        di as error "cpsmatch: internal error (matching engine)"
        exit 197
    }
    if `__plugin' quietly gen double _weight = .
    else {
        quietly gen double _weight = _treated if _support == 1
        if "`ate'" != "" quietly replace _weight = 0 if _treated == 1 & _support == 1
    }
    char _weight[Type] "aweight"
    label variable _weight "cpsmatch: weight of matched controls"
    foreach v of local outcome {
        if `__plugin' quietly gen double _`v' = .
        else if "`ate'" == "" quietly gen double _`v' = 0 if _support == 1 & _treated == 1
        else quietly gen double _`v' = 0 if _support == 1
        label variable _`v' "cpsmatch: value of `v' of match(es)"
    }

    if "`warnings'" != "" & "`metric'" == "pscore" {
        _cpsmatch_warn, method(`method') `work_ate'
    }

    * Neighbor identifiers: psmatch2's _id is the position in its matching
    * order (support, _treated, score, observation); _n# hold the positions
    * of each row's matches.
    if `nearest' {
        quietly gen long _id = .
        label variable _id "cpsmatch: Identifier (ID)"
        forvalues n = 1/`neighbor' {
            quietly gen long _n`n' = .
            label variable _n`n' "cpsmatch: ID of nearest neighbor nr. `n'"
        }
        local n1 _n1
        local nk _n`neighbor'
    }
    * _nn: matches per row; for radius and kernel matching on the score,
    * the controls with a positive weight
    if `nearest' | `__plugin' {
        if `__plugin' quietly gen long _nn = .
        else quietly gen long _nn = 0 if _support == 1
        label variable _nn "cpsmatch: # matched neighbors"
    }
    if `nearest' & "`metric'" == "pscore" {
        quietly gen double _pdif = .
        label variable _pdif "cpsmatch: abs(pscore - pscore[nearest neighbor])"
    }
    if `ai' > 0 {
        foreach v of local outcome {
            capture drop _self_`v'
            quietly gen double _self_`v' = .
            label variable _self_`v' "cpsmatch: matched value of `v' (T-T & C-C)"
        }
    }

    if `__do_timing' {
        quietly _ctools_clock __ct_now
        local __ct_elapsed91 = `__ct_elapsed91' + (`__ct_now' - `__ct_start91')
        quietly _ctools_clock __ct_start92
    }

    * =========================================================================
    * MATCH
    * =========================================================================

    * Sufficient statistics come from the plugin for its ATT-only results;
    * ate, ai() and the smoothed and Mata paths compute psmatch2's table here.
    local plugin_stats = ("`engine'" == "plugin" & "`ate'" == "" & `ai' == 0 & `smooth' == 0 & !`radius_all')
    if "`engine'" == "plugin" {
        local threads_code ""
        if `threads' > 0 local threads_code "threads(`threads')"
        local method_code = cond(`nearest', 0, cond("`radius'" != "", 1, 2))
        local with_replace = ("`noreplacement'" == "")
        local handle_ties = ("`ties'" != "")
        local desc_order = ("`descending'" != "")
        * The support before matching is _support itself, which the plugin
        * reads before it stores; without common or trim every row of the
        * sample is on support and the plugin needs no such variable.
        local __presup = cond("`common'" != "" | `trim' < 100, "presup=3", "")

        * A plugin call replaces _support, so later calls (one per outcome)
        * get the support before matching back from a copy. presupport
        * marks the rows on support before matching.
        local presupport `touse'
        if "`__presup'" != "" {
            tempvar presupport
            quietly gen byte `presupport' = _support
        }
        capture scalar drop _cpsmatch_ai_pass

        * Several outcomes share one match (it depends on the score alone):
        * the plugin matches once per outcome and stores its matched values.
        local outcomes `outcome'
        if "`outcomes'" == "" local outcomes _cpsmatch_no_outcome_
        local k_outcome = 0
        foreach cur_outcome of local outcomes {
            if "`cur_outcome'" == "_cpsmatch_no_outcome_" local cur_outcome
            local ++k_outcome
            local cur_match : word `k_outcome' of `matchvars'
            * Plugin variables: treatment _pscore _support _weight _nn
            * [_id _n1.._nK _pdif] [outcome _<outcome> [_self_<outcome>]]
            * _treated
            local allvars "`depvar' _pscore _support _weight _nn"
            local flags "`__presup' nnvar=5"
            local id_idx = 0
            local p = 5
            if `nearest' {
                local id_idx = 6
                local allvars "`allvars' _id `n1'-`nk' _pdif"
                local flags "`flags' ids=7 pdif=`=7 + `neighbor''"
                local p = 7 + `neighbor'
            }
            local outcome_idx = 0
            if "`cur_outcome'" != "" {
                local outcome_idx = `p' + 1
                local allvars "`allvars' `cur_match' _`cur_outcome'"
                local flags "`flags' matched=`=`p' + 2'"
                local p = `p' + 2
                if `ai' > 0 {
                    local allvars "`allvars' _self_`cur_outcome'"
                    local flags "`flags' self=`=`p' + 1' ai=`ai'"
                    local p = `p' + 1
                    if "`altvariance'" != "" local flags "`flags' altvar"
                }
            }
            * (without an outcome, the support counts and missing moments)
            if `plugin_stats' local flags "`flags' stats"
            if "`work_ate'" != "" local flags "`flags' ate"
            * The arms are the stored _treated (a code it cannot hold is
            * missing there); the treatment itself marks the sample rows.
            local allvars "`allvars' _treated"
            local flags "`flags' arm=`=`p' + 1'"
            * Format: treat_idx pscore_idx outcome_idx weight_idx id_idx
            *         support_idx nobs method neighbor caliper common kernel
            *         bwidth with_replace ties desc, then flags. The sample
            *         marker goes last, in place of if `touse'.
            local args "1 2 `outcome_idx' 4 `id_idx' 3"
            local args "`args' `nobs' `method_code' `neighbor' `caliper' 0"
            local args "`args' `kernelcode' `bwidth' `with_replace' `handle_ties' `desc_order'"
            if `k_outcome' > 1 & "`__presup'" != "" quietly replace _support = `presupport'
            plugin call ctools_plugin `allvars' `touse', "cpsmatch `threads_code' `args' `flags' touse"
        }
        * ai(#) above 1 with noreplacement: only a row's first same-arm match
        * is taken without replacement, the others are the rows around it
        * (psmatch2). A plugin that says so sets _cpsmatch_ai_pass to 2; an
        * older one takes them all without replacement and runs out of rows.
        if `ai' > 1 & "`noreplacement'" != "" & "`outcome'" != "" {
            capture confirm scalar _cpsmatch_ai_pass
            if _rc {
                di as error "cpsmatch: rebuild the ctools plugin alongside the ado files"
                exit 601
            }
        }
        * The correction also needs, for each column of the score model, the
        * mean over a row's ai() same-arm matches: the same pass as
        * _self_<outcome>, with the column in the outcome's place. The other
        * outputs of these calls repeat the ones above and go to scratch
        * variables.
        local psselfx
        if `pscorr' & "`outcome'" != "" {
            tempvar xsupport xweight xnn xmatched
            quietly gen byte `xsupport' = .
            quietly gen double `xweight' = .
            quietly gen long `xnn' = .
            quietly gen double `xmatched' = .
            foreach x of local psx {
                tempvar xself
                quietly gen double `xself' = .
                local psselfx `psselfx' `xself'
                local allvars "`depvar' _pscore `xsupport' `xweight' `xnn' `x' `xmatched' `xself' _treated"
                local flags "nnvar=5 matched=7 self=8 ai=`ai' ate arm=9"
                local args "1 2 6 4 0 3 `nobs' `method_code' `neighbor' `caliper' 0"
                local args "`args' `kernelcode' `bwidth' `with_replace' `handle_ties' `desc_order'"
                plugin call ctools_plugin `allvars' `touse', "cpsmatch `threads_code' `args' `flags' touse"
            }
            drop `xsupport' `xweight' `xnn' `xmatched'
        }
        * Rows outside the sample follow it in psmatch2's order: by score
        * (they have one when it was given with pscore()), then as they are
        if `nearest' & `nobs' < _N {
            mata: _cpsmatch_outside_ids("`touse'", `nobs', "`descending'" != "")
        }
        * psmatch2's same-arm pass on an arm with one row on support before
        * matching reads the row after it in the matching order, whatever
        * that row is (the first row of the other arm, a row with another
        * treatment code, off support or outside the sample), and takes it as
        * the row's match if it has a score within the caliper. cpsmatch
        * does the same. When no row follows, psmatch2 stops with a Mata
        * error (r(3301)); cpsmatch stops here.
        if `nearest' & `ai' > 0 {
            * the outcomes and, for the correction, the score model's columns
            local self_from `outcome'
            local self_to
            foreach v of local outcome {
                local self_to `self_to' _self_`v'
            }
            if "`psselfx'" != "" {
                local self_from `self_from' `psx'
                local self_to `self_to' `psselfx'
            }
            local self_cal = cond(`caliper' == 0, ., `caliper')
            forvalues arm = 0/1 {
                quietly count if `presupport' == 1 & _treated == `arm'
                if r(N) != 1 continue
                mata: _cpsmatch_self_single("`presupport'", `arm', `self_cal', ///
                    "`altvariance'" != "", "`self_from'", "`self_to'")
                if !`self_found' {
                    di as error "cpsmatch: ai() has no same-arm match for the only " ///
                        cond(`arm', "treated", "untreated") " observation on support"
                    exit 2001
                }
            }
        }
    }
    else {
        _cpsmatch_mata, method(`method') metric(`metric') neighbor(`neighbor') ///
            caliper(`caliper') kerneltype(`kerneltype') bwidth(`bwidth') ai(`ai') ///
            outcome(`outcome') mahalanobis(`mahalanobis') w(`w') order(`order') ///
            touse(`touse') `ate' `altvariance'
    }
    if `nearest' {
        * psmatch2 compresses _id and gives _n# its type
        quietly compress _id
        local idtype : type _id
        if "`idtype'" != "long" {
            forvalues n = 1/`neighbor' {
                quietly recast `idtype' _n`n'
            }
        }
    }

    * psmatch2's radius without a caliper is uniform-kernel matching with a
    * missing bandwidth: abs(dif) <= . holds for every observation of the
    * data, so each target spreads its weight equally over all _N of them
    * (the donors' mean is still their plain mean). On support, every
    * observation's _weight thus gains (number of targets) / _N.
    if `radius_all' {
        quietly count if _support == 1 & (_treated == 1 | "`ate'" != "")
        local n_targets = r(N)
        quietly replace _weight = cond("`ate'" != "", 0, _treated) + `n_targets' / _N if _support == 1
        quietly replace _weight = . if _weight == 0 | _support == 0
    }

    if `__do_timing' {
        quietly _ctools_clock __ct_now
        local __ct_elapsed92 = `__ct_elapsed92' + (`__ct_now' - `__ct_start92')
        quietly _ctools_clock __ct_start94
    }

    * =========================================================================
    * TREATMENT EFFECTS AND TABLE
    * =========================================================================

    local shown_header = 0
    local k_outcome = 0
    * r(table), as psmatch2 builds it: a column per outcome and effect
    tempname rtable rcolumn psq
    local rtable_eq
    local rtable_names
    local rtable_crit = invnormal(1 - (100 - c(level))/200)
    foreach cur_outcome of local outcome {
        local ++k_outcome
        local att = .
        local att_se = .
        local att_t = .
        local atu = .
        local ate_v = .
        local seatu = .
        local seate = .
        * The plugin's statistics are those of the last outcome's call (the
        * earlier outcomes' are replaced), and stand for psmatch2's sums when
        * every control of the sample is on support, as in an ATT-only match.
        local from_plugin = (`plugin_stats' & `k_outcome' == `: word count `outcome'')
        if `from_plugin' local from_plugin = (scalar(_cpsmatch_n0_off) == 0)
        if `from_plugin' {
            local y1_bar = scalar(_cpsmatch_y1)
            local n1 = scalar(_cpsmatch_n1)
            local var1 = scalar(_cpsmatch_var1)
            local y0_bar = scalar(_cpsmatch_y0)
            local n0 = scalar(_cpsmatch_n0)
            local var1_u = scalar(_cpsmatch_var1_u)
            local n1_u = scalar(_cpsmatch_n1_u)
            local var0_u = scalar(_cpsmatch_var0_u)
            local n0_u = scalar(_cpsmatch_n0_u)
            if `n1' > 0 & `n0' > 0 {
                local att = `y1_bar' - `y0_bar'
                * psmatch2's default ATT variance: the unweighted outcome
                * variances of the treated and of the controls on support
                * (here all the controls), the latter times the sum of the
                * controls' squared weights.
                local sum_weight_squared = scalar(_cpsmatch_w2)
                local att_se = sqrt(`var1'/`n1' + `var0_u'*`sum_weight_squared'/`n1'^2)
            }
            local y1_unmatched = scalar(_cpsmatch_y1_u)
            local y0_unmatched = scalar(_cpsmatch_y0_u)
            * psmatch2 displays the OLS SE from regressing outcome on
            * treatment: its pooled residual variance, where a
            * singleton group contributes zero SSE.
            local sse1_u = cond(`n1_u' > 1, (`n1_u'-1)*`var1_u', 0)
            local sse0_u = cond(`n0_u' > 1, (`n0_u'-1)*`var0_u', 0)
            local pooled_var = (`sse1_u' + `sse0_u')/(`n1_u' + `n0_u' - 2)
            local unmatched_se = sqrt(`pooled_var'*(1/`n1_u' + 1/`n0_u'))
        }
        else {
            local skipse = ""
            if `smooth' == 2 | ("`method'" == "llr" & "`kerneltype'" == "tricube") local skipse skipse
            _cpsmatch_effects `cur_outcome', neighbor(`neighbor') ai(`ai') `work_ate' ///
                `samplevar' `altvariance' `skipse'
            local y1_bar = r(m1t)
            local y0_bar = r(m0t)
            local att = r(att)
            local att_se = r(seatt)
            local y1_unmatched = r(u1u)
            local y0_unmatched = r(u0u)
            local unmatched_se = r(seols)
            if "`work_ate'" != "" {
                local m1u = r(m1u)
                local m0u = r(m0u)
                local atu = r(atu)
                local ate_v = r(ate)
                local seatu = r(seatu)
                local seate = r(seate)
                local n1_on = r(n1)
                local n0_on = r(n0)
            }
            if `pscorr' {
                * Abadie and Imbens (2016): the variances less the part the
                * estimated score takes away, plus (ATT and ATU) the part
                * the effect's dependence on the score adds
                local fixed_att = `att_se'
                local fixed_atu = `seatu'
                local fixed_ate = `seate'
                mata: st_matrix("`psq'", _cpsmatch_ai2016("`cur_outcome'", "`psx'", "`psselfx'", ///
                    "`psdens'", "`psV'", "`pskeep'", `pscons' > 0, `ai', `neighbor', ///
                    `n0_on', `n1_on', `att', `atu'))
                local v_ate = `fixed_ate'^2 - `psq'[1, 1]
                local v_att = `fixed_att'^2 - `psq'[1, 2] + `psq'[1, 3]
                local v_atu = `fixed_atu'^2 - `psq'[1, 4] + `psq'[1, 5]
                foreach e in ate att atu {
                    if `v_`e'' >= -1e-10 & `v_`e'' < 0 local v_`e' = 0
                    if `v_`e'' >= 0 local v_`e' = sqrt(`v_`e'')
                    else {
                        di as text "Warning: corrected `=upper("`e'")' variance < -1e-10 for `cur_outcome'; SE set to missing."
                        local v_`e' = .
                    }
                }
                local att_se = `v_att'
                local seatu = `v_atu'
                local seate = `v_ate'
            }
        }
        if !missing(`att_se') & `att_se' > 0 local att_t = `att' / `att_se'
        local effects att
        if "`ate'" != "" local effects att atu ate
        foreach e of local effects {
            local b = cond("`e'" == "att", `att', cond("`e'" == "atu", `atu', `ate_v'))
            local s = cond("`e'" == "att", `att_se', cond("`e'" == "atu", `seatu', `seate'))
            matrix `rcolumn' = J(9, 1, .)
            matrix `rcolumn'[1, 1] = `b'
            matrix `rcolumn'[2, 1] = `s'
            if `s' < . & `s' > 0 {
                matrix `rcolumn'[3, 1] = (`b')/`s'
                matrix `rcolumn'[4, 1] = 2*normal(-abs((`b')/`s'))
                matrix `rcolumn'[5, 1] = `b' - `rtable_crit'*`s'
                matrix `rcolumn'[6, 1] = `b' + `rtable_crit'*`s'
            }
            matrix `rcolumn'[8, 1] = `rtable_crit'
            matrix `rcolumn'[9, 1] = 0
            matrix `rtable' = nullmat(`rtable'), `rcolumn'
            local rtable_eq `rtable_eq' `cur_outcome'
            local rtable_names `rtable_names' `=upper("`e'")'
        }

        * ATT results table (psmatch2 format)
        if `att' != . | "`ate'" != "" {
            local unmatched_diff = `y1_unmatched' - `y0_unmatched'
            local unmatched_t = `unmatched_diff' / `unmatched_se'
            if !`shown_header' {
                di as text "{hline 78}"
                di as text "        Variable     Sample |    Treated     Controls   Difference         S.E.   T-stat"
                di as text "{hline 28}+{hline 59}"
                local shown_header = 1
            }
            di as text %16s abbrev("`cur_outcome'", 16) "  Unmatched |" as result %11.6g `y1_unmatched' %12.6g `y0_unmatched' %12.6g `unmatched_diff' %12.6g `unmatched_se' %9.2f `unmatched_t'
            di as text "                        ATT |" as result %11.6g `y1_bar' %12.6g `y0_bar' %12.6g `att' %12.6g `att_se' %9.2f `att_t'
            if "`ate'" != "" {
                di as text "                        ATU |" as result %11.6g `m0u' %12.6g `m1u' %12.6g `atu' %12.6g `seatu' %9.2f `atu'/`seatu'
                di as text "                        ATE |" as result _col(53) %12.6g `ate_v' %12.6g `seate' %9.2f `ate_v'/`seate'
            }
            di as text "{hline 28}+{hline 59}"
        }
        * psmatch2's per-outcome results, where the names fit in 32 characters
        foreach r in att seatt atu ate seatu seate {
            if "`ate'" == "" & inlist("`r'", "atu", "ate", "seatu", "seate") continue
            local value = cond("`r'" == "att", `att', cond("`r'" == "seatt", `att_se', ///
                cond("`r'" == "atu", `atu', cond("`r'" == "ate", `ate_v', ///
                cond("`r'" == "seatu", `seatu', `seate')))))
            if strlen("`r'_`cur_outcome'") <= 32 return scalar `r'_`cur_outcome' = `value'
            local last_`r' = `value'
        }
        * debug: the standard errors before the correction and its parts
        if `pscorr' & "`debug'" != "" {
            local parts "seatt_ai_fixed qTminus qTplus"
            if "`ate'" != "" local parts "`parts' seate_ai_fixed seatu_ai_fixed qA qUminus qUplus"
            foreach r of local parts {
                if "`r'" == "seatt_ai_fixed" local value = `fixed_att'
                else if "`r'" == "seate_ai_fixed" local value = `fixed_ate'
                else if "`r'" == "seatu_ai_fixed" local value = `fixed_atu'
                else if "`r'" == "qA" local value = `psq'[1, 1]
                else if "`r'" == "qTminus" local value = `psq'[1, 2]
                else if "`r'" == "qTplus" local value = `psq'[1, 3]
                else if "`r'" == "qUminus" local value = `psq'[1, 4]
                else local value = `psq'[1, 5]
                if strlen("`r'_`cur_outcome'") <= 32 return scalar `r'_`cur_outcome' = `value'
            }
        }
        local last_att_t = `att_t'
    }
    if "`outcome'" != "" {
        matrix rownames `rtable' = b se z pvalue ll ul df crit eform
        matrix coleq `rtable' = `rtable_eq'
        matrix colnames `rtable' = `rtable_names'
        return matrix table = `rtable'
    }
    if `shown_header' {
        if `ai' == 0 {
            if `last_seatt' != . {
                di as text "Note: S.E. does not take into account that the propensity score is estimated."
                di as text "      Use ai(#) for Abadie-Imbens analytical matching standard errors."
            }
        }
        else if `pscorr' di as text "Note: Population AI S.E. adjusted for estimated propensity scores."
        else {
            if "`samplevar'" != "" di as text "Note: Sample AI S.E."
            else di as text "Note: Population AI S.E."
            if `psfixnote' di as text "Note: AI S.E. does not take into account that the propensity score is estimated because the AI(2016) correction is not available for this specification."
        }
    }

    * psmatch2 matched the controls to the treated for the correction alone:
    * without ate it clears their matched outcomes, neighbors and _nn again
    * and gives the treated their weight of 1 (_pdif stays for the controls).
    if `pscorr' & "`ate'" == "" {
        * (psmatch2 loops over the outcomes here, and stops without any)
        if "`outcome'" == "" {
            di as error "varlist required"
            exit 100
        }
        foreach v of local outcome {
            quietly replace _`v' = . if _treated == 0
        }
        forvalues n = 1/`neighbor' {
            quietly replace _n`n' = . if _treated == 0
        }
        quietly replace _nn = . if _treated == 0
        quietly replace _weight = 1 if _treated == 1 & _support == 1
        quietly replace _weight = . if _support == 0 | _weight == 0
    }

    if `__do_timing' {
        quietly _ctools_clock __ct_now
        local __ct_elapsed94 = `__ct_elapsed94' + (`__ct_now' - `__ct_start94')
    }

    * Treatment assignment table (psmatch2 format); the plugin counts its
    * rows, the Mata paths are counted here
    di as text ""
    if "`engine'" == "plugin" {
        local n_untreated_on = scalar(_cpsmatch_n_controls_on)
        local n_untreated_off = scalar(_cpsmatch_n_controls_off)
        local n_treated_on = scalar(_cpsmatch_n_matched)
        local n_treated_off = scalar(_cpsmatch_n_off_support)
    }
    else {
        quietly count if _treated == 0 & _support == 1
        local n_untreated_on = r(N)
        quietly count if _treated == 0 & _support == 0
        local n_untreated_off = r(N)
        quietly count if _treated == 1 & _support == 1
        local n_treated_on = r(N)
        quietly count if _treated == 1 & _support == 0
        local n_treated_off = r(N)
    }
    local n_untreated_total = `n_untreated_on' + `n_untreated_off'
    local n_treated_total = `n_treated_on' + `n_treated_off'
    local n_on_total = `n_untreated_on' + `n_treated_on'
    local n_off_total = `n_untreated_off' + `n_treated_off'
    local n_total = `n_on_total' + `n_off_total'

    di as text "           | cpsmatch:"
    di as text " cpsmatch: |   Common"
    di as text " Treatment |  support"
    if `n_off_total' > 0 {
        di as text "assignment | Off suppo  On suppor |     Total"
        di as text "{hline 11}+{hline 22}+{hline 10}"
        di as text " Untreated |" as result %10.0f `n_untreated_off' %10.0f `n_untreated_on' as text " |" as result %10.0f `n_untreated_total'
        di as text "   Treated |" as result %10.0f `n_treated_off' %10.0f `n_treated_on' as text " |" as result %10.0f `n_treated_total'
        di as text "{hline 11}+{hline 22}+{hline 10}"
        di as text "     Total |" as result %10.0f `n_off_total' %10.0f `n_on_total' as text " |" as result %10.0f `n_total'
    }
    else {
        di as text "assignment | On suppor |     Total"
        di as text "{hline 11}+{hline 11}+{hline 10}"
        di as text " Untreated |" as result %10.0f `n_untreated_on' as text " |" as result %10.0f `n_untreated_total'
        di as text "   Treated |" as result %10.0f `n_treated_on' as text " |" as result %10.0f `n_treated_total'
        di as text "{hline 11}+{hline 11}+{hline 10}"
        di as text "     Total |" as result %10.0f `n_on_total' as text " |" as result %10.0f `n_total'
    }

    di as text ""

    * =========================================================================
    * TIMING OUTPUT (verbose)
    * =========================================================================

    if `__do_timing' {
        quietly _ctools_clock __ct_now
        local __ct_elapsed90 = `__ct_elapsed90' + (`__ct_now' - `__ct_start90')
        local __time_total = `__ct_elapsed90'
        local __time_preplugin = `__ct_elapsed91'
        local __time_plugin = `__ct_elapsed92'
        local __time_att = `__ct_elapsed94'

        di as text "{hline 55}"
        di as text "cpsmatch timing breakdown:"
        di as text "{hline 55}"
        if "`engine'" == "plugin" {
            di as text "  C plugin internals:"
            di as text "    Load data:              " as result %8.4f _cpsmatch_time_load " sec"
            di as text "    Group separation:       " as result %8.4f _cpsmatch_time_setup " sec"
            di as text "    Sort/index build:       " as result %8.4f _cpsmatch_time_sort " sec"
            di as text "    Matching:               " as result %8.4f _cpsmatch_time_match " sec"
            di as text "    Outcome statistics:     " as result %8.4f _cpsmatch_time_stats " sec"
            di as text "    Store results:          " as result %8.4f _cpsmatch_time_store " sec"
            di as text "  {hline 53}"
            di as text "    C plugin total:         " as result %8.4f _cpsmatch_time_total " sec"
            di as text ""
        }
        else {
            di as text "    Matching (Mata):        " as result %8.4f `__time_plugin' " sec"
        }
        di as text "  Stata overhead:"
        di as text "    Pre-plugin setup:       " as result %8.4f `__time_preplugin' " sec"
        di as text "    Post-plugin (ATT):      " as result %8.4f `__time_att' " sec"
        di as text "  {hline 53}"
        local __stata_overhead = `__time_preplugin' + `__time_att'
        di as text "    Stata overhead total:   " as result %8.4f `__stata_overhead' " sec"
        di as text "{hline 55}"
        di as text "    Wall clock total:       " as result %8.4f `__time_total' " sec"
        di as text "{hline 55}"

        capture local __threads_max = _cpsmatch_threads_max
        if _rc == 0 & "`engine'" == "plugin" {
            capture local __openmp_enabled = _cpsmatch_openmp_enabled
            if _rc != 0 local __openmp_enabled = 0
            di as text ""
            di as text "  Thread diagnostics:"
            di as text "    OpenMP enabled:         " as result %8.0f `__openmp_enabled'
            di as text "    Max threads available:  " as result %8.0f `__threads_max'
            di as text "{hline 55}"
        }
    }

    * =========================================================================
    * RETURN RESULTS
    * =========================================================================

    local n_treated = `n_treated_on' + `n_treated_off'
    local n_controls = `n_untreated_on' + `n_untreated_off'
    local cmin = .
    local cmax = .
    if "`engine'" != "plugin" {
        capture confirm variable _pscore
        if !_rc {
            quietly summarize _pscore if _treated == 0, meanonly
            if r(N) {
                local cmin = r(min)
                local cmax = r(max)
            }
        }
    }
    return scalar N = `nobs'
    return scalar n_treated = `n_treated'
    return scalar n_controls = `n_controls'
    return scalar n_matched = `n_treated_on'
    return scalar n_off_support = `n_treated_off'
    if "`engine'" == "plugin" {
        return scalar common_min = scalar(_cpsmatch_common_min)
        return scalar common_max = scalar(_cpsmatch_common_max)
    }
    else {
        return scalar common_min = `cmin'
        return scalar common_max = `cmax'
    }

    * psmatch2 returns the effects of the last outcome, missing when
    * undefined (each outcome's are returned above)
    if "`outcome'" != "" {
        return scalar att = `last_att'
        return scalar att_se = `last_seatt'
        return scalar seatt = `last_seatt'
        return scalar att_t = `last_att_t'
        if "`ate'" != "" {
            foreach r in atu ate seatu seate {
                return scalar `r' = `last_`r''
            }
        }
    }

    * psmatch2's record of the ai() standard errors
    return scalar ai = `ai_given'
    if `ai_given' > 0 {
        return local ai_vartype = cond("`samplevar'" != "", "sample", "population")
        return scalar ai_pscorr = `pscorr'
        if `pscorr' return local ai_pscorr_method "AI2016_ps_decomposition"
    }

    local method_name = cond("`radius'" != "", "radius", cond("`method'" == "neighbor", "nearest", "`method'"))
    if `smooth' == 1 local method_name "llr"
    if `smooth' == 2 local method_name "spline"
    return local method "`method_name'"
    return local metric "`metric'"
    return local treatvar "`depvar'"
    if "`outcome'" != "" {
        return local outcome "`outcome'"
        * psmatch2's names for the outcome list and the matching covariates
        return local depvar "`outcome'"
    }
    if "`treatvar'" != "" return local exog "`treatvar'"
end


* psmatch2's common support: the treated outside the controls' score range,
* or with trim(#) the treated where the controls' score density is zero or
* in its lowest # percent (with ate, the controls likewise). The score is the
* one estimated or given, before odds.
program define _cpsmatch_support
    syntax varname, level(real) [ate]
    if `level' == 100 {
        summarize `varlist' if _treated == 0, meanonly
        replace _support = 0 if (`varlist' < r(min) | `varlist' > r(max)) & _treated == 1
        if "`ate'" != "" {
            summarize `varlist' if _treated == 1, meanonly
            replace _support = 0 if (`varlist' < r(min) | `varlist' > r(max)) & _treated == 0
        }
    }
    else {
        _cpsmatch_trim `varlist', level(`level') arm(1)
        if "`ate'" != "" _cpsmatch_trim `varlist', level(`level') arm(0)
    }
end

program define _cpsmatch_trim
    syntax varname, level(real) arm(integer)
    tempvar at density
    kdensity `varlist' if _treated == (1 - `arm'), nograph at(`varlist') generate(`at' `density')
    replace _support = 0 if `density' == 0 & _treated == `arm'
    if `level' > 0 {
        _pctile `density' if _treated == `arm', percentiles(`level')
        replace _support = 0 if `density' < r(r1) & _treated == `arm'
    }
end


* psmatch2's warnings: controls (with ate, also treated) with equal scores
program define _cpsmatch_warn, sortpreserve
    syntax, method(string) [ate]
    sort _treated _pscore
    capture by _treated _pscore: assert _N == 1 if _treated == 0 & _support == 1
    if !_rc & "`ate'" != "" {
        capture by _treated _pscore: assert _N == 1 if _treated == 1 & _support == 1
    }
    if _rc & "`method'" == "neighbor" {
        di as result "There are observations with identical propensity score values."
        di as result "The sort order of the data could affect your results."
        di as result "Make sure that the sort order is random before calling cpsmatch."
    }
end


* Treatment effects of one outcome from the generated variables, with
* psmatch2's formulas: the ATT (and ATU, ATE), their approximate standard
* errors, or with ai(#) the Abadie-Imbens standard errors.
program define _cpsmatch_effects, rclass
    syntax varname, [neighbor(integer 1) ai(integer 0) ate samplevar altvariance skipse]
    local v `varlist'
    tempname att atu effect seatt seatu seate m1t m0t m1u m0u
    quietly {
        replace _`v' = . if _support == 0
        capture replace _self_`v' = . if _support == 0

        summarize `v' if _treated == 1
        local u1u = r(mean)
        local n1u = r(N)
        local var1u = r(Var)
        summarize `v' if _treated == 0
        local u0u = r(mean)
        local n0u = r(N)
        local var0u = r(Var)

        summarize `v' if _treated == 1 & _support == 1, meanonly
        scalar `m1t' = r(mean)
        local n1 = r(N)
        summarize _`v' if _treated == 1 & _support == 1, meanonly
        scalar `m0t' = r(mean)
        scalar `att' = `m1t' - `m0t'
        scalar `atu' = .
        scalar `effect' = .
        scalar `m1u' = .
        scalar `m0u' = .
        local n0 = 0
        if "`ate'" != "" {
            summarize _`v' if _treated == 0 & _support == 1, meanonly
            scalar `m1u' = r(mean)
            summarize `v' if _treated == 0 & _support == 1, meanonly
            scalar `m0u' = r(mean)
            local n0 = r(N)
            scalar `atu' = `m1u' - `m0u'
            scalar `effect' = `att'*`n1'/(`n0' + `n1') + `atu'*`n0'/(`n0' + `n1')
        }

        scalar `seatu' = .
        scalar `seate' = .
        if `ai' != 0 {
            * Abadie and Imbens (2006, eq. 14), or Abadie et al. (2004, p. 303)
            * with altvariance: the variance of the population effect
            * (their Theorem 7), or with samplevar that of the sample effect
            * (Theorem 6). As psmatch2, the terms are generated in the
            * default storage type.
            tempvar w shat vt vu ve
            generate `w' = max(_weight, 0)
            if "`altvariance'" == "" generate `shat' = (`ai'/(`ai' + 1))*(`v' - _self_`v')^2
            else generate `shat' = _self_`v'
            if "`samplevar'" != "" {
                generate `vt' = `shat'*(_treated - (1 - _treated)*`w')^2
                if "`ate'" != "" {
                    generate `vu' = `shat'*((1 - _treated) - _treated*`w')^2
                    generate `ve' = `shat'*(1 + `w')^2
                }
            }
            else {
                generate `vt' = max(0, _treated*(`v' - _`v' - `att')^2) ///
                    + (1 - _treated)*(`w'^2 - `w'/`neighbor')*`shat'
                if "`ate'" != "" {
                    generate `vu' = max(0, (1 - _treated)*(_`v' - `v' - `atu')^2) ///
                        + _treated*(`w'^2 - `w'/`neighbor')*`shat'
                    generate `ve' = max(0, (_treated*(`v' - _`v') + (1 - _treated)*(_`v' - `v') ///
                        - `effect')^2) + (`w'^2 + 2*`w' - `w'/`neighbor')*`shat'
                }
            }
            summarize `vt' if _support == 1, meanonly
            scalar `seatt' = sqrt(r(sum))/`n1'
            if "`ate'" != "" {
                summarize `vu' if _support == 1, meanonly
                scalar `seatu' = sqrt(r(sum))/`n0'
                summarize `ve' if _support == 1, meanonly
                scalar `seate' = sqrt(r(sum))/(`n0' + `n1')
            }
        }
        else {
            * psmatch2's approximate standard errors: the outcome variances
            * of the treated and of the controls on support, each arm's
            * times the sum of its squared weights (for the ATE, of one plus
            * the weight, over the observations used as a match: an unused
            * one has a missing weight).
            tempvar w2
            summarize `v' if _treated == 1 & _support == 1
            local var1 = r(Var)
            summarize `v' if _treated == 0 & _support == 1
            local var0 = r(Var)
            generate double `w2' = _weight^2 if _support == 1
            summarize `w2' if _treated == 0, meanonly
            local wtot = r(sum)
            if "`skipse'" != "" scalar `seatt' = .
            else scalar `seatt' = sqrt(`var1'/`n1' + `var0'*`wtot'/`n1'^2)
            if "`ate'" != "" & "`skipse'" == "" {
                summarize `w2' if _treated == 1, meanonly
                scalar `seatu' = sqrt(`var0'/`n0' + `var1'*r(sum)/`n0'^2)
                replace `w2' = (1 + _weight)^2 if _support == 1
                summarize `w2' if _treated == 1, meanonly
                local wtot = r(sum)
                summarize `w2' if _treated == 0, meanonly
                scalar `seate' = sqrt((`var1'*`wtot' + `var0'*r(sum))/(`n0' + `n1')^2)
            }
        }
    }
    * The Unmatched row: the OLS standard error of the outcome on treatment
    local sse1 = cond(`n1u' > 1, (`n1u' - 1)*`var1u', 0)
    local sse0 = cond(`n0u' > 1, (`n0u' - 1)*`var0u', 0)
    return scalar seols = sqrt((`sse1' + `sse0')/(`n1u' + `n0u' - 2)*(1/`n1u' + 1/`n0u'))
    return scalar u1u = `u1u'
    return scalar u0u = `u0u'
    return scalar m1t = `m1t'
    return scalar m0t = `m0t'
    return scalar m1u = `m1u'
    return scalar m0u = `m0u'
    return scalar att = `att'
    return scalar atu = `atu'
    return scalar ate = `effect'
    return scalar seatt = `seatt'
    return scalar seatu = `seatu'
    return scalar seate = `seate'
    return scalar n1 = `n1'
    return scalar n0 = `n0'
end


* Mahalanobis-metric matching (nearest neighbor, kernel, llr) and local
* linear regression matching, over psmatch2's matching order.
program define _cpsmatch_mata
    syntax, method(string) metric(string) neighbor(integer) caliper(real) ///
        order(varname) touse(varname) [kerneltype(string) bwidth(real 0.06) ai(integer 0) ///
        outcome(varlist) mahalanobis(string) w(string) ate altvariance]

    tempvar msup
    quietly gen byte `msup' = -_support
    sort `msup' _treated `order'
    local nearest = ("`method'" == "neighbor")
    if `nearest' {
        quietly replace _id = _n
    }
    * Blocks of the order: the supported controls, then the supported treated
    quietly count if _support == 1 & _treated < 0
    local c0 = r(N) + 1
    quietly count if _treated == 0 & _support == 1
    local c1 = `c0' + r(N) - 1
    quietly count if _support == 1 & _treated > 0 & _treated < 1
    local t0 = `c1' + r(N) + 1
    quietly count if _treated == 1 & _support == 1
    local t1 = `t0' + r(N) - 1
    local cal = cond(`caliper' == 0, ., `caliper')

    local xvars
    tempname W
    if "`metric'" == "mahalanobis" {
        fvrevar `mahalanobis'
        local originals `r(varlist)'
        tempname scales
        matrix `scales' = J(1, `: word count `originals'', 1)
        local j = 0
        foreach v of local originals {
            local ++j
            tempvar centered
            quietly generate double `centered' = `v'
            mata: _cpsmatch_center("`centered'", "`touse'", "`scales'", `j')
            local xvars `xvars' `centered'
        }
        quietly _rmdcoll _treated `xvars' if `touse', expand
        fvrevar `r(varlist)'
        local xvars `r(varlist)'
        if "`w'" == "" {
            tempname XX0 XX1
            quietly matrix accum `XX0' = `xvars' if _treated == 0, deviations noconstant
            quietly matrix accum `XX1' = `xvars' if _treated == 1, deviations noconstant
            quietly count if _treated <= 1 & _support == 1
            matrix `W' = syminv((`XX0' + `XX1')/(r(N) - 2))
        }
        else {
            matrix `W' = diag(`scales') * `w' * diag(`scales')
        }
    }
    else {
        local xvars _pscore
        matrix `W' = (1)
    }
    local matched
    foreach v of local outcome {
        local matched `matched' _`v'
    }
    local is_mahal = ("`metric'" == "mahalanobis")

    if `nearest' {
        local ids
        forvalues n = 1/`neighbor' {
            local ids `ids' _n`n'
        }
        if "`ate'" != "" {
            mata: _cpsmatch_nn(`c0', `c1', `t0', `t1', `neighbor', `cal', "`xvars'", "`W'", ///
                "`ids'", "`outcome'", "`matched'", 0, 1, 0)
            quietly replace _support = 0 if _n1 >= . & _treated == 0
        }
        mata: _cpsmatch_nn(`t0', `t1', `c0', `c1', `neighbor', `cal', "`xvars'", "`W'", ///
            "`ids'", "`outcome'", "`matched'", 0, "`ate'" != "", 0)
        quietly replace _support = 0 if _n1 >= . & _treated == 1
        if `ai' > 0 {
            * An arm of one row has no same-arm match, and psmatch2 stops
            * with a Mata error (r(3301)) unless it skips the row (a control
            * that no treated row was matched to, without ate).
            foreach b in c t {
                if ``b'0' == ``b'1' & ("`ate'" != "" | _weight[``b'0'] != 0) {
                    di as error "cpsmatch: ai() has no same-arm match for the only " ///
                        cond("`b'" == "t", "treated", "untreated") " observation on support"
                    exit 2001
                }
            }
            local self
            foreach v of local outcome {
                quietly replace _self_`v' = 0 if _support == 1
                local self `self' _self_`v'
            }
            local alt = ("`altvariance'" != "")
            mata: _cpsmatch_nn(`c0', `c1', `c0', `c1', `ai', `cal', "`xvars'", "`W'", ///
                "`ids'", "`outcome'", "`self'", 1, "`ate'" != "", `alt')
            mata: _cpsmatch_nn(`t0', `t1', `t0', `t1', `ai', `cal', "`xvars'", "`W'", ///
                "`ids'", "`outcome'", "`self'", 1, "`ate'" != "", `alt')
        }
    }
    else {
        local llr = ("`method'" == "llr")
        mata: _cpsmatch_smooth(`c0', `c1', `t0', `t1', `llr', `is_mahal', "`xvars'", "`W'", ///
            "`kerneltype'", `bwidth', "`outcome'", "`matched'", "`ate'" != "")
    }
    quietly replace _weight = . if _weight == 0 | _support == 0
end


mata:

// Center before rank screening; retain units for a user-supplied metric.
void _cpsmatch_center(string scalar name, string scalar sample,
    string scalar scales, real scalar j)
{
    real colvector x, keep, z
    real matrix s
    real scalar origin, spread
    x = st_data(., name)
    keep = selectindex(st_data(., sample) :!= 0)
    origin = x[keep[1]]
    z = x :- origin
    spread = max(abs(z[keep]))
    if (spread >= .) _error(498)
    if (spread == 0) spread = 1
    st_store(., name, z :/ spread)
    s = st_matrix(scales)
    s[j] = spread
    st_matrix(scales, s)
}

// The k smallest distances of a (distance, position) matrix, by Hoare's
// selection (as psmatch2, whose list order this keeps); rows beyond k are
// not ordered.
real matrix _cpsmatch_select(real matrix a, real scalar k)
{
    real scalar lo, hi, i, j
    real scalar pivot
    real rowvector swap

    lo = 1
    hi = rows(a)
    while (lo < hi) {
        pivot = a[k, 1]
        i = lo
        j = hi
        do {
            while (a[i, 1] < pivot) i++
            while (pivot < a[j, 1]) j--
            if (i <= j) {
                swap = a[i, .]
                a[i, .] = a[j, .]
                a[j, .] = swap
                i++
                j--
            }
        } while (i <= j)
        if (j < k) lo = i
        if (k < i) hi = j
    }
    return(a[|1, 1 \ k, 2|])
}

// Nearest-neighbor matching on the quadratic metric x'Wx of xvars: targets
// t0..t1 to donors c0..c1 (positions of the sorted data). A match needs a
// squared distance within caliper. Between arms (self = 0) it records _n#,
// _nn, _weight (+1/neighbors per match) and the matched outcomes; with
// self = 1 (ai) the targets skip themselves, controls unused without ate
// are skipped, and the outcomes go to the _self_ variables, as conditional
// variances with altvar.
void _cpsmatch_nn(real scalar t0, real scalar t1, real scalar c0, real scalar c1,
    real scalar neighbors, real scalar caliper, string scalar xvars, string scalar wname,
    string scalar idvars, string scalar outvars, string scalar matchvars,
    real scalar self, real scalar ate, real scalar altvar)
{
    real matrix X, W, IDS, Y, M, d, best, delta
    real colvector q, weight, nn, support, pick
    real rowvector mean
    real scalar i, m, n, k, nout, first, last

    if (t1 < t0 | c1 < c0) {
        if (!self & t1 >= t0) {
            st_view(support, ., "_support")
            support[|t0 \ t1|] = J(t1 - t0 + 1, 1, 0)
        }
        return
    }
    st_view(X, ., tokens(xvars))
    W = st_matrix(wname)
    st_view(weight, ., "_weight")
    st_view(support, ., "_support")
    if (!self) {
        st_view(IDS, ., tokens(idvars))
        st_view(nn, ., "_nn")
    }
    nout = cols(tokens(outvars))
    if (nout) {
        st_view(Y, ., tokens(outvars))
        st_view(M, ., tokens(matchvars))
    }
    first = min((t0, c0))
    last = max((t1, c1))
    for (i = t0; i <= t1; i++) {
        if (self & !ate & weight[i] == 0) continue
        delta = X[|c0, 1 \ c1, .|] :- X[i, .]
        d = rowsum((delta * W) :* delta)
        d = d, (c0::c1)
        if (self) d = select(d, d[., 2] :!= i)
        k = min((neighbors, rows(d)))
        if (k < 1) pick = J(0, 1, .)
        else {
            best = _cpsmatch_select(d, k)
            pick = select(best[., 2], best[., 1] :<= caliper)
        }
        n = rows(pick)
        if (!self) nn[i] = n
        for (m = 1; m <= n; m++) {
            if (!self) {
                IDS[i, m] = pick[m]
                weight[pick[m]] = weight[pick[m]] + 1 / neighbors
            }
            if (nout) M[i, .] = M[i, .] + Y[pick[m], .] :/ n
        }
        if (altvar & self & nout) {
            mean = (n :* M[i, .] + Y[i, .]) :/ (n + 1)
            M[i, .] = (Y[i, .] - mean):^2
            for (m = 1; m <= n; m++) M[i, .] = M[i, .] + (Y[pick[m], .] - mean):^2
            M[i, .] = M[i, .] :/ n
        }
        if (n < 1 & !self) support[i] = 0
    }
}

// psmatch2's _id outside the sample: after the sample's rows (first of
// them), by score (by its negative with descending), then in the data's
// order; rows without a score come last
void _cpsmatch_outside_ids(string scalar touse, real scalar first, real scalar descending)
{
    real colvector at, key

    at = selectindex(st_data(., touse) :== 0)
    if (!rows(at)) return
    key = st_data(at, "_pscore")
    if (descending) key = -key
    st_store(at[order((key, at), (1, 2))], "_id", first :+ (1::rows(at)))
}

// ai(#) on an arm with one row on support before matching (presup): the
// row's same-arm match is the row after it in the matching order (_id), if
// that row has a score within caliper, as psmatch2's pass has it. The row's
// dstvars then get that row's srcvars (with altvar the conditional variance
// of Abadie et al. 2004, p. 303, of the two rows), when the row is still on
// support. The local self_found tells whether a row follows at all.
void _cpsmatch_self_single(string scalar presup, real scalar arm, real scalar caliper,
    real scalar altvar, string scalar srcvars, string scalar dstvars)
{
    real colvector id, row, next
    real rowvector own, value, center
    string rowvector src

    row = selectindex(st_data(., presup) :== 1 :& st_data(., "_treated") :== arm)
    id = st_data(., "_id")
    next = selectindex(id :== id[row[1]] + 1)
    st_local("self_found", strofreal(rows(next) > 0))
    src = tokens(srcvars)
    if (!rows(next) | !cols(src)) return
    if (st_data(row[1], "_support") != 1) return
    if (!(abs(st_data(row[1], "_pscore") - st_data(next[1], "_pscore")) < caliper)) return
    value = st_data(next[1], src)
    if (altvar) {
        own = st_data(row[1], src)
        center = (value + own) :/ 2
        value = (own - center):^2 + (value - center):^2
    }
    st_store(row[1], tokens(dstvars), value)
}

// psmatch2's kernel weights of distances d with bandwidth bw, normalized to
// sum to one; missing outside the kernel's support
real colvector _cpsmatch_kernel(real colvector d, string scalar kernel, real scalar bw)
{
    real colvector u, w, inside

    u = d :/ bw
    inside = abs(d) :<= bw
    if (kernel == "normal") w = normalden(u)
    else if (kernel == "epan") w = 1 :- u:^2
    else if (kernel == "biweight") w = (1 :- u:^2):^2
    else if (kernel == "uniform") w = J(rows(d), 1, 1)
    else w = (1 :- abs(u):^3):^3
    if (kernel != "normal") w = w :/ inside
    return(w :/ quadsum(w))
}

// Kernel and local linear regression matching: each target (the supported
// treated, with ate first the supported controls) against the supported
// rows of the other arm, on the score difference or the quadratic metric.
void _cpsmatch_smooth(real scalar c0, real scalar c1, real scalar t0, real scalar t1,
    real scalar llr, real scalar mahal, string scalar xvars, string scalar wname,
    string scalar kernel, real scalar bw, string scalar outvars, string scalar matchvars,
    real scalar ate)
{
    real matrix X, W, Y, M, delta
    real colvector q, weight, support, treated, donors, d, w, ok, y, dd, ww, keep
    real scalar i, j, n, s, v, mean, sw, b, nout, first, arm, cons, sxx, sxy, xbar, ybar

    st_view(X, ., tokens(xvars))
    W = st_matrix(wname)
    st_view(weight, ., "_weight")
    st_view(support, ., "_support")
    st_view(treated, ., "_treated")
    nout = cols(tokens(outvars))
    if (nout) {
        st_view(Y, ., tokens(outvars))
        st_view(M, ., tokens(matchvars))
    }
    first = ate ? c0 : t0
    for (i = first; i <= t1; i++) {
        if (i > c1 & i < t0) continue
        arm = (i <= c1)
        donors = selectindex(support :== 1 :& treated :== arm)
        if (rows(donors) == 0) {
            // as psmatch2, llr leaves support to its regressions
            if (!llr | nout) support[i] = 0
            continue
        }
        if (mahal) {
            delta = X[donors, .] :- X[i, .]
            d = rowsum((delta * W) :* delta)
        }
        else d = X[donors, 1] :- X[i, 1]
        w = _cpsmatch_kernel(d, kernel, bw)
        ok = selectindex(w :< .)
        if (!llr) {
            if (rows(ok) == 0) {
                support[i] = 0
                continue
            }
            for (j = 1; j <= rows(ok); j++) {
                weight[donors[ok[j]]] = weight[donors[ok[j]]] + w[ok[j]]
            }
            for (v = 1; v <= nout; v++) {
                y = Y[donors[ok], v]
                ww = w[ok]
                ww = select(ww, y :< .)
                y = select(y, y :< .)
                if (rows(y)) {
                    mean = quadsum(ww :* y) / quadsum(ww)
                    if (mean < .) M[i, v] = mean
                }
            }
            continue
        }
        // llr: the local linear weights (psmatch2's, from summarize with
        // aweights, which drops zero weights) and the intercept of a
        // weighted regression on d
        dd = select(d[ok], w[ok] :> 0)
        ww = select(w[ok], w[ok] :> 0)
        n = rows(ww)
        sw = quadsum(ww)
        s = quadsum(ww :* dd)
        mean = s / sw
        b = n > 1 ? quadsum(ww :* (dd :- mean):^2) / sw : .
        for (j = 1; j <= rows(ok); j++) {
            weight[donors[ok[j]]] = weight[donors[ok[j]]] + w[ok[j]] * (b + s^2 - s * d[ok[j]]) / b
        }
        for (v = 1; v <= nout; v++) {
            y = Y[donors[ok], v]
            keep = (y :< .) :& (w[ok] :> 0)
            ww = select(w[ok], keep)
            dd = select(d[ok], keep)
            y = select(y, keep)
            if (rows(y) < 2) {
                support[i] = 0
                continue
            }
            sw = quadsum(ww)
            xbar = quadsum(ww :* dd) / sw
            ybar = quadsum(ww :* y) / sw
            sxx = quadsum(ww :* (dd :- xbar):^2)
            sxy = quadsum(ww :* (dd :- xbar) :* (y :- ybar))
            cons = sxx > 0 ? ybar - (sxy / sxx) * xbar : ybar
            M[i, v] = cons
        }
    }
}

// Abadie and Imbens (2016): what estimating the propensity score changes in
// the variances of nearest-neighbor matching on it, as psmatch2 4.1.0
// computes it. With c a row vector over the score model's coefficients and V
// their variance, each returned term is c V c': (ATE, ATT minus, ATT plus,
// ATU minus, ATU plus). The ATE variance loses its term; the ATT and ATU
// variances lose the "minus" term and gain the "plus" term, the effect's
// derivative with respect to the coefficients.
//
// Per observation on support: X is its row of the model's columns (the
// constant last), f the derivative of its score p with respect to the index,
// and C1, C0 estimates of the covariance of X with the outcome under
// treatment and under control at its score. The covariance in the
// observation's own arm is ai/(ai + 1) (X - mean X)(Y - mean Y) over its ai()
// same-arm matches; in the other arm it is the average of that over its
// matches there. mu1 - mu0 is the same-arm matched mean outcome against the
// other arm's.
real rowvector _cpsmatch_ai2016(string scalar yname, string scalar xnames,
    string scalar selfnames, string scalar densname, string scalar vname,
    string scalar keep, real scalar hascons, real scalar ai, real scalar neighbors,
    real scalar n0, real scalar n1, real scalar att, real scalar atu)
{
    real matrix X, own, other, c1, c0, V
    real colvector sel, y, sy, my, f, p, q, treated, id, rowof, count, donor, pick, pos, ok
    real colvector delta, aT, aU
    real rowvector idx, cA, cT, dT, cU, dU
    real scalar n, m

    idx = strtoreal(tokens(keep))
    V = st_matrix(vname)[idx, idx]
    sel = selectindex(st_data(., "_support") :== 1)
    n = rows(sel)
    if (n == 0) return(J(1, 5, .))
    p = st_data(sel, "_pscore")
    if (any(p :<= 0) | any(p :>= 1)) {
        _error(3498, "cpsmatch: propensity score outside (0,1) in the ai() correction")
    }
    q = 1 :- p
    y = st_data(sel, yname)
    sy = st_data(sel, "_self_" + yname)
    my = st_data(sel, "_" + yname)
    f = st_data(sel, densname)
    treated = st_data(sel, "_treated") :== 1
    X = st_data(sel, tokens(xnames))
    own = (ai / (ai + 1)) :* (X - st_data(sel, tokens(selfnames))) :* (y - sy)

    // the other arm's covariance: the average over the matches, by their ids
    id = st_data(sel, "_id")
    rowof = J(st_nobs(), 1, .)
    rowof[id] = 1::n
    other = J(n, cols(X), 0)
    count = J(n, 1, 0)
    for (m = 1; m <= neighbors; m++) {
        donor = st_data(sel, "_n" + strofreal(m))
        pick = selectindex(donor :< .)
        if (rows(pick) == 0) continue
        count[pick] = count[pick] :+ 1
        pos = rowof[donor[pick]]
        ok = selectindex(pos :< .)
        if (rows(ok)) other[pick[ok], .] = other[pick[ok], .] + own[pos[ok], .]
    }
    if (any(count :== 0)) {
        _error(3498, "cpsmatch: an observation on support has no match in the ai() correction")
    }
    other = other :/ count

    c1 = other + treated :* (own - other)
    c0 = own + treated :* (other - own)
    delta = (2 :* treated :- 1) :* (sy - my)
    aT = delta :- att
    aU = delta :- atu
    if (hascons) {
        X = X, J(n, 1, 1)
        c1 = c1, J(n, 1, 0)
        c0 = c0, J(n, 1, 0)
    }
    cA = colsum(f :* (c1 :/ p + c0 :/ q)) :/ (n0 + n1)
    cT = colsum(f :* (X :* aT + c1 + (p :/ q) :* c0)) :/ n1
    dT = colsum(f :* (X :* aT + c1 - c0)) :/ n1
    cU = colsum(f :* (X :* aU - c0 - (q :/ p) :* c1)) :/ n0
    dU = colsum(f :* (X :* aU + c1 - c0)) :/ n0
    return((cA * V * cA', cT * V * cT', dT * V * dT', cU * V * cU', dU * V * dU'))
}

end
