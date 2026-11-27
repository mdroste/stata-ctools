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

    * The plugin makes _treated and _support itself (common support
    * included) when nothing reads them before matching: no smoothing, no
    * trim() and no Mata engine. They are then created without filling, as
    * the plugin writes every observation; paths that stop with an error
    * before matching keep generating them here.
    local __fast = !`__mata' & !`willsmooth' & `trim' == 100
    if ("`common'" != "" | `trim' < 100) & "`score'" == "" local __fast = 0
    if `__fast' {
        quietly count if `touse'
        local nobs = r(N)
        if `nobs' == 0 local __fast = 0
    }
    * (fast: _treated and _support are made with the outputs below, in one
    * step, as each addition of variables rewrites the dataset)
    capture label drop _treated
    label define _treated 0 "Untreated" 1 "Treated"
    if !`__fast' {
        quietly gen byte _treated = `depvar' if `touse'
        label variable _treated "cpsmatch: Treatment assignment"
        label values _treated _treated
    }

    * =========================================================================
    * COMMON SUPPORT
    * =========================================================================

    if ("`common'" != "" | `trim' < 100) & "`score'" == "" {
        di as error "With option 'common' a propensity score is needed. Provide one"
        di as error "with option 'pscore()' or estimate one. See the help file for more details."
        exit 198
    }
    capture label drop _support
    label define _support 0 "Off support" 1 "On support"
    if !`__fast' {
        quietly gen byte _support = 1 if `touse'
        label variable _support "cpsmatch: Common support"
        label values _support _support
    }
    if ("`common'" != "" | `trim' < 100) & !`__fast' {
        if !inrange(`trim', 0, 100) {
            di as error "Trim level out of range"
            exit 198
        }
        quietly _cpsmatch_support `score', level(`trim') `ate'
    }

    if !`__fast' {
        quietly count if `touse'
        local nobs = r(N)
    }
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
    * The plugin writes every observation of its outputs, so they are
    * created in one step without filling, in psmatch2's order and with
    * their final types: _id and _n# get the type compress gives _id, whose
    * values are 1.._N.
    if `__plugin' {
        if `ai' > 0 {
            foreach v of local outcome {
                capture drop _self_`v'
            }
        }
        local __idtype = cond(_N <= 100, "byte", cond(_N <= 32740, "int", "long"))
        local __types double
        local __names _weight
        if `__fast' {
            local __types byte byte `__types'
            local __names _treated _support `__names'
        }
        foreach v of local outcome {
            local __types `__types' double
            local __names `__names' _`v'
        }
        if `nearest' {
            local __types `__types' `__idtype'
            local __names `__names' _id
            forvalues n = 1/`neighbor' {
                local __types `__types' `__idtype'
                local __names `__names' _n`n'
            }
        }
        local __types `__types' long
        local __names `__names' _nn
        if `nearest' {
            local __types `__types' double
            local __names `__names' _pdif
        }
        if `ai' > 0 {
            foreach v of local outcome {
                local __types `__types' double
                local __names `__names' _self_`v'
            }
        }
        * the same-arm means of the score model's columns, for the ai()
        * correction (temporary, after the outputs)
        local psselfx
        if `pscorr' & "`outcome'" != "" & "`psx'" != "" {
            foreach x of local psx {
                tempvar xself
                local psselfx `psselfx' `xself'
                local __types `__types' double
            }
        }
        mata: (void) st_addvar(tokens("`__types'"), tokens("`__names' `psselfx'"), 1)
        if `nearest' & "`__idtype'" != "long" format _id _n1-_n`neighbor' %12.0g
        if `__fast' {
            label variable _treated "cpsmatch: Treatment assignment"
            label values _treated _treated
            label variable _support "cpsmatch: Common support"
            label values _support _support
        }
    }
    else {
        quietly gen double _weight = _treated if _support == 1
        if "`ate'" != "" quietly replace _weight = 0 if _treated == 1 & _support == 1
    }
    char _weight[Type] "aweight"
    label variable _weight "cpsmatch: weight of matched controls"
    foreach v of local outcome {
        if !`__plugin' {
            if "`ate'" == "" quietly gen double _`v' = 0 if _support == 1 & _treated == 1
            else quietly gen double _`v' = 0 if _support == 1
        }
        label variable _`v' "cpsmatch: value of `v' of match(es)"
    }

    * (the plugin checks the scores for ties itself)
    if "`warnings'" != "" & "`metric'" == "pscore" & !`__plugin' {
        _cpsmatch_warn, method(`method') `work_ate'
    }

    * Neighbor identifiers: psmatch2's _id is the position in its matching
    * order (support, _treated, score, observation); _n# hold the positions
    * of each row's matches.
    if `nearest' {
        if !`__plugin' quietly gen long _id = .
        label variable _id "cpsmatch: Identifier (ID)"
        forvalues n = 1/`neighbor' {
            if !`__plugin' quietly gen long _n`n' = .
            label variable _n`n' "cpsmatch: ID of nearest neighbor nr. `n'"
        }
        local n1 _n1
        local nk _n`neighbor'
    }
    * _nn: matches per row; for radius and kernel matching on the score,
    * the controls with a positive weight
    if `nearest' | `__plugin' {
        if !`__plugin' quietly gen long _nn = 0 if _support == 1
        label variable _nn "cpsmatch: # matched neighbors"
    }
    if `nearest' & "`metric'" == "pscore" {
        if !`__plugin' quietly gen double _pdif = .
        label variable _pdif "cpsmatch: abs(pscore - pscore[nearest neighbor])"
    }
    if `ai' > 0 {
        foreach v of local outcome {
            if !`__plugin' {
                capture drop _self_`v'
                quietly gen double _self_`v' = .
            }
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
        * The support before matching: with trim() or smoothing, _support as
        * Stata made it (the plugin reads it before it stores); otherwise the
        * plugin makes it, imposing common support itself (with ate also on
        * the controls).
        local __presup = cond(!`__fast' & ("`common'" != "" | `trim' < 100), "presup=3", "")
        local __common = 0
        if `__fast' & "`common'" != "" local __common = cond("`ate'" != "", 2, 1)
        capture scalar drop _cpsmatch_ai_pass

        * One call matches once for every outcome (the match depends on the
        * score alone) and, for the ai() correction, for the score model's
        * columns. Plugin variables: treatment _pscore _support _weight _nn
        * [_id _n1.._nK _pdif] [the outcomes matched, _<outcome> [_self_<outcome>]
        * [the outcomes, when smoothed ones are matched]] [the score model's
        * columns and their same-arm means] [the support score] _treated,
        * then the sample marker. The plugin writes every observation (fill)
        * and _id outside the sample (outids), and with ai() an arm of one row
        * on support takes the row after it in the matching order, whatever
        * that row is, as psmatch2's pass does (single; with no row after it,
        * psmatch2 stops with a Mata error, r(3301), and cpsmatch with r(2001)).
        local ny : word count `outcome'
        local allvars "`depvar' _pscore _support _weight _nn"
        local flags "`__presup' nnvar=5 fill"
        local id_idx = 0
        local p = 5
        if `nearest' {
            local id_idx = 6
            local allvars "`allvars' _id `n1'-`nk' _pdif"
            local flags "`flags' ids=7 pdif=`=7 + `neighbor'' outids"
            local p = 7 + `neighbor'
        }
        local outcome_idx = 0
        if `ny' > 0 {
            local outcome_idx = `p' + 1
            local __mlist
            local __slist
            foreach v of local outcome {
                local __mlist `__mlist' _`v'
                local __slist `__slist' _self_`v'
            }
            local allvars "`allvars' `matchvars' `__mlist'"
            local flags "`flags' ny=`ny' matched=`=`p' + `ny' + 1'"
            local p = `p' + 2*`ny'
            if `ai' > 0 {
                local allvars "`allvars' `__slist'"
                local flags "`flags' self=`=`p' + 1'"
                local p = `p' + `ny'
                if "`altvariance'" != "" local flags "`flags' altvar"
            }
            if `smooth' {
                local allvars "`allvars' `outcome'"
                local flags "`flags' yeff=`=`p' + 1'"
                local p = `p' + `ny'
            }
            * what _cpsmatch_effects summarizes, for each outcome
            local flags "`flags' effects"
            if "`samplevar'" != "" local flags "`flags' samplevar"
        }
        if `ai' > 0 local flags "`flags' ai=`ai' single"
        * The correction also needs, for each column of the score model, the
        * mean over a row's ai() same-arm matches: the same pass as
        * _self_<outcome>, with the column in the outcome's place.
        if "`psselfx'" != "" {
            local __nx : word count `psx'
            local allvars "`allvars' `psx' `psselfx'"
            local flags "`flags' nx=`__nx' xcol=`=`p' + 1' xself=`=`p' + `__nx' + 1'"
            local p = `p' + 2*`__nx'
        }
        if `__common' & "`score'" != "_pscore" {
            local allvars "`allvars' `score'"
            local flags "`flags' supscore=`=`p' + 1'"
            local p = `p' + 1
        }
        * (without an outcome, the support counts and missing moments)
        if `plugin_stats' local flags "`flags' stats"
        if "`work_ate'" != "" local flags "`flags' ate"
        if "`warnings'" != "" local flags "`flags' warn"
        if `radius_all' local flags "`flags' radall"
        * The arms are the stored _treated (a code it cannot hold is
        * missing there; mkarm: the plugin makes it); the treatment itself
        * marks the sample rows.
        local allvars "`allvars' _treated"
        local flags "`flags' arm=`=`p' + 1'"
        if `__fast' local flags "`flags' mkarm"
        * Format: treat_idx pscore_idx outcome_idx weight_idx id_idx
        *         support_idx nobs method neighbor caliper common kernel
        *         bwidth with_replace ties desc, then flags. The sample
        *         marker goes last, in place of if `touse'.
        local args "1 2 `outcome_idx' 4 `id_idx' 3"
        local args "`args' `nobs' `method_code' `neighbor' `caliper' `__common'"
        local args "`args' `kernelcode' `bwidth' `with_replace' `handle_ties' `desc_order'"
        capture noisily plugin call ctools_plugin `allvars' `touse', "cpsmatch `threads_code' `args' `flags' touse"
        local __rc = _rc
        if `__rc' & `__rc' != 2001 {
            * the state the generated variables had before the plugin's stores
            foreach v of local __names {
                quietly replace `v' = .
            }
            if `__fast' {
                quietly replace _treated = cond(`touse', `depvar', .)
                quietly replace _support = cond(`touse', 1, .)
            }
            exit `__rc'
        }
        if `__rc' {
            * (r(2001): psmatch2's stop comes before _id is compressed)
            if `nearest' & "`__idtype'" != "long" quietly recast long _id _n1-_n`neighbor'
            exit `__rc'
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
        * psmatch2's warnings: controls (with ate, also treated) with equal scores
        if "`warnings'" != "" & "`method'" == "neighbor" {
            if scalar(_cpsmatch_dup) {
                di as result "There are observations with identical propensity score values."
                di as result "The sort order of the data could affect your results."
                di as result "Make sure that the sort order is random before calling cpsmatch."
            }
        }
    }
    else {
        _cpsmatch_mata, method(`method') metric(`metric') neighbor(`neighbor') ///
            caliper(`caliper') kerneltype(`kerneltype') bwidth(`bwidth') ai(`ai') ///
            outcome(`outcome') mahalanobis(`mahalanobis') w(`w') order(`order') ///
            touse(`touse') `ate' `altvariance'
    }
    if `nearest' & !`__plugin' {
        * psmatch2 compresses _id and gives _n# its type (the plugin's
        * outputs are made with that type)
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
    * (the plugin does this for its own matching)
    if `radius_all' & !`__plugin' {
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
            * the plugin summarized this outcome's matching results
            local __from = cond(`__plugin', "from(`k_outcome')", "")
            _cpsmatch_effects `cur_outcome', neighbor(`neighbor') ai(`ai') `work_ate' ///
                `samplevar' `altvariance' `skipse' `__from'
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
                * The plugin computes the terms over the observations on
                * support; without ate it then clears the controls' matched
                * outcome (and after the last outcome _n#, _nn and _weight),
                * as psmatch2 leaves them when it matched the controls for
                * the correction alone.
                local __aiflags
                if "`ate'" == "" {
                    local __aiflags clear
                    if `k_outcome' == `: word count `outcome'' local __aiflags clear final
                }
                plugin call ctools_plugin _pscore `psdens' `cur_outcome' _self_`cur_outcome' ///
                    _`cur_outcome' _treated _id `n1'-`nk' `psx' `psselfx' _support _weight _nn ///
                    `touse', "cpsmatch ai2016 `: word count `psx'' `neighbor' `ai' `=`pscons' > 0' `n0_on' `n1_on' `att' `atu' `psV' `pskeep' `__aiflags'"
                matrix `psq' = (scalar(_cpsmatch_psq1), scalar(_cpsmatch_psq2), ///
                    scalar(_cpsmatch_psq3), scalar(_cpsmatch_psq4), scalar(_cpsmatch_psq5))
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
    * (the plugin's correction calls did that above, outcome by outcome)
    if `pscorr' & "`ate'" == "" {
        * (psmatch2 loops over the outcomes here, and stops without any)
        if "`outcome'" == "" {
            di as error "varlist required"
            exit 100
        }
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
    syntax varname, [neighbor(integer 1) ai(integer 0) ate samplevar altvariance skipse ///
        from(integer 0)]
    local v `varlist'
    tempname att atu effect seatt seatu seate m1t m0t m1u m0u
    * from(#): the plugin has summarized outcome # of its call. Its scalars
    * _cpsmatch_e#_* are what the summarize and generate steps below give
    * (Stata's sums are exact sums rounded once, which the plugin
    * reproduces), and _cpsmatch_e_* the weight sums, which no outcome
    * changes; the formulas on them stay as they are, locals included.
    local E _cpsmatch_e`from'_
    quietly {
        if `from' {
            local u1u = scalar(`E'u1u)
            local n1u = scalar(`E'n1u)
            local var1u = scalar(`E'var1u)
            local u0u = scalar(`E'u0u)
            local n0u = scalar(`E'n0u)
            local var0u = scalar(`E'var0u)
            scalar `m1t' = scalar(`E'm1t)
            local n1 = scalar(`E'n1)
            scalar `m0t' = scalar(`E'm0t)
        }
        else {
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
        }
        scalar `att' = `m1t' - `m0t'
        scalar `atu' = .
        scalar `effect' = .
        scalar `m1u' = .
        scalar `m0u' = .
        local n0 = 0
        if "`ate'" != "" {
            if `from' {
                scalar `m1u' = scalar(`E'm1u)
                scalar `m0u' = scalar(`E'm0u)
                local n0 = scalar(`E'n0)
            }
            else {
                summarize _`v' if _treated == 0 & _support == 1, meanonly
                scalar `m1u' = r(mean)
                summarize `v' if _treated == 0 & _support == 1, meanonly
                scalar `m0u' = r(mean)
                local n0 = r(N)
            }
            scalar `atu' = `m1u' - `m0u'
            scalar `effect' = `att'*`n1'/(`n0' + `n1') + `atu'*`n0'/(`n0' + `n1')
        }

        scalar `seatu' = .
        scalar `seate' = .
        if `ai' != 0 & `from' {
            * (the plugin's sums of the terms generated below)
            scalar `seatt' = sqrt(scalar(`E'sumvt))/`n1'
            if "`ate'" != "" {
                scalar `seatu' = sqrt(scalar(`E'sumvu))/`n0'
                scalar `seate' = sqrt(scalar(`E'sumve))/(`n0' + `n1')
            }
        }
        else if `ai' != 0 {
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
            if `from' {
                local var1 = scalar(`E'var1)
                local var0 = scalar(`E'var0)
                local wtot = scalar(_cpsmatch_e_w2c)
            }
            else {
                summarize `v' if _treated == 1 & _support == 1
                local var1 = r(Var)
                summarize `v' if _treated == 0 & _support == 1
                local var0 = r(Var)
                generate double `w2' = _weight^2 if _support == 1
                summarize `w2' if _treated == 0, meanonly
                local wtot = r(sum)
            }
            if "`skipse'" != "" scalar `seatt' = .
            else scalar `seatt' = sqrt(`var1'/`n1' + `var0'*`wtot'/`n1'^2)
            if "`ate'" != "" & "`skipse'" == "" & `from' {
                scalar `seatu' = sqrt(`var0'/`n0' + `var1'*scalar(_cpsmatch_e_w2t)/`n0'^2)
                local wtot = scalar(_cpsmatch_e_w1t)
                scalar `seate' = sqrt((`var1'*`wtot' + `var0'*scalar(_cpsmatch_e_w1c))/(`n0' + `n1')^2)
            }
            else if "`ate'" != "" & "`skipse'" == "" {
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
        * The plugin's matching loops, operation for operation as psmatch2's
        * Mata code: targets, donors, neighbors, caliper, metric, self, ate,
        * altvariance, then the counts of metric columns and outcomes
        local nx : word count `xvars'
        local nout : word count `outcome'
        if "`ate'" != "" {
            plugin call ctools_plugin `xvars' _weight _support _nn `ids' `outcome' `matched', ///
                "cpsmatch mnn `c0' `c1' `t0' `t1' `neighbor' `cal' `W' 0 1 0 `nx' `nout'"
            quietly replace _support = 0 if _n1 >= . & _treated == 0
        }
        plugin call ctools_plugin `xvars' _weight _support _nn `ids' `outcome' `matched', ///
            "cpsmatch mnn `t0' `t1' `c0' `c1' `neighbor' `cal' `W' 0 `=("`ate'" != "")' 0 `nx' `nout'"
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
            local ate1 = ("`ate'" != "")
            plugin call ctools_plugin `xvars' _weight _support `outcome' `self', ///
                "cpsmatch mnn `c0' `c1' `c0' `c1' `ai' `cal' `W' 1 `ate1' `alt' `nx' `nout'"
            plugin call ctools_plugin `xvars' _weight _support `outcome' `self', ///
                "cpsmatch mnn `t0' `t1' `t0' `t1' `ai' `cal' `W' 1 `ate1' `alt' `nx' `nout'"
        }
    }
    else {
        local llr = ("`method'" == "llr")
        local kcode = 0
        if "`kerneltype'" == "normal" local kcode = 1
        else if "`kerneltype'" == "biweight" local kcode = 2
        else if "`kerneltype'" == "uniform" local kcode = 3
        else if "`kerneltype'" == "tricube" local kcode = 4
        local nx : word count `xvars'
        local nout : word count `outcome'
        plugin call ctools_plugin `xvars' _weight _support _treated `outcome' `matched', ///
            "cpsmatch msmooth `c0' `c1' `t0' `t1' `llr' `is_mahal' `W' `kcode' `bwidth' `=("`ate'" != "")' `nx' `nout'"
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

end

* ctools-helpers:begin (generated by scripts/sync_release.py from src/ado; edit those files)
* Private monotonic timestamp in seconds; no public Stata timer IDs.
program define _ctools_clock
    version 14.1
    args target
    tempname held
    _return hold `held'
    capture noisily {
        * Once registered, this plugin stays loaded for the Stata session.
        * Resolve its path and check its protocol only on first use (or
        * after discard), rather than searching the ado-path at every lap.
        capture plugin call ctools_clock_plugin, "clock"
        if _rc {
            _ctools_load
            capture program ctools_clock_plugin, plugin using("`__ctools_plugin'")
            if _rc != 0 & _rc != 110 exit _rc
            plugin call ctools_clock_plugin, "clock"
        }
    }
    local rc = _rc
    _return restore `held'
    if `rc' exit `rc'
    c_local `target' `ctools_clock_value'
end

* Shared plugin loader and identity contract
program define _ctools_load
    version 14.1
    * Not rclass, and findfile's r() is undone, so the caller's r() survives;
    * identity comes back through c_local.
    local __os = c(os)
    local __machine = c(machine_type)
    local __is_mac = 0
    if "`__os'" == "MacOSX" {
        local __is_mac = 1
    }
    else if strpos(lower("`__machine'"), "mac") > 0 {
        local __is_mac = 1
    }

    local __plugin = ""
    if "`__os'" == "Windows" {
        local __plugin "ctools_windows.plugin"
    }
    else if `__is_mac' {
        local __is_arm = 0
        if strpos(lower("`__machine'"), "apple") > 0 | strpos(lower("`__machine'"), "arm") > 0 | strpos(lower("`__machine'"), "silicon") > 0 {
            local __is_arm = 1
        }
        if `__is_arm' == 0 & !strpos(lower("`__machine'"), "intel") & !strpos(lower("`__machine'"), "x86") {
            tempfile __archfile
            quietly shell uname -m > "`__archfile'" 2>&1
            tempname __fh
            file open `__fh' using "`__archfile'", read text
            file read `__fh' __archline
            file close `__fh'
            capture erase "`__archfile'"
            if strpos("`__archline'", "arm64") > 0 {
                local __is_arm = 1
            }
        }
        if `__is_arm' {
            local __plugin "ctools_mac_arm.plugin"
        }
        else {
            local __plugin "ctools_mac_x86.plugin"
        }
    }
    else if "`__os'" == "Unix" {
        local __plugin "ctools_linux.plugin"
    }
    else {
        local __plugin "ctools.plugin"
    }

    * The ado/plugin argument protocol, stamped from config/release.json.
    local __ado_protocol "4"
    local __reinstall `"net install ctools, from("https://raw.githubusercontent.com/mdroste/stata-ctools/main/build") replace"'

    tempname __rhold
    _return hold `__rhold'
    capture findfile `__plugin'
    local __found = (_rc == 0)
    local __path `"`r(fn)'"'
    _return restore `__rhold'
    local __wanted "`__plugin'"
    if !`__found' & "`__plugin'" != "ctools.plugin" local __plugin "ctools.plugin"
    capture program ctools_plugin, plugin using("`__plugin'")
    local __load_rc = _rc
    if `__load_rc' != 0 & `__load_rc' != 110 {
        di as error "ctools: could not load the plugin `__plugin'"
        di as error `"  Stata reports c(os) "`__os'" and c(machine_type) "`__machine'"; ctools selected `__wanted'."'
        if `__found' di as error `"  The plugin file is `__path'."'
        else di as error "  `__wanted' was not found on the ado-path; reinstall ctools."
        di as error "  Stata's plugin loader reported:"
        capture noisily program ctools_plugin, plugin using("`__plugin'")
        if "`__os'" == "Unix" & !`__is_mac' {
            di as error "  The distribution Linux plugin needs x86-64 and glibc 2.17 or later; OpenMP is linked statically."
            di as error "  Run ldd -r on the plugin to inspect missing libraries or symbols."
        }
        else if `__is_mac' {
            di as error "  On macOS, check that the plugin is the notarized ctools release and is not quarantined."
        }
        di as error "  See https://github.com/mdroste/stata-ctools/blob/main/docs/INSTALLATION.md#troubleshooting"
        exit `__load_rc'
    }
    local ctools_plugin_hdfe_screen
    local ctools_plugin_hdfe_markout
    local ctools_plugin_hdfe_fvlevels
    capture plugin call ctools_plugin, "version"
    local __id_rc = _rc
    if `__id_rc' | "`ctools_plugin_protocol'" != "`__ado_protocol'" {
        local __has = cond("`ctools_plugin_protocol'" == "", "none", "`ctools_plugin_protocol'")
        di as error "ctools: the plugin `__plugin' does not match the installed ado files."
        local __ident = cond("`ctools_plugin_version'" == "", "", " (plugin version `ctools_plugin_version', revision `ctools_plugin_revision')")
        di as error "  The ado files need plugin protocol `__ado_protocol'; the loaded plugin reports `__has'`__ident'."
        di as error "  Reinstall ctools so its ado files and plugins come from the same release:"
        di as error `"    `__reinstall'"'
        di as error "  Then restart Stata, because a plugin stays loaded until Stata exits."
        exit 9
    }
    c_local __ctools_plugin "`__plugin'"
    c_local __ctools_plugin_version "`ctools_plugin_version'"
    c_local __ctools_plugin_revision "`ctools_plugin_revision'"
    c_local __ctools_plugin_protocol "`ctools_plugin_protocol'"
    c_local __ctools_plugin_hdfe_screen = ("`ctools_plugin_hdfe_screen'" == "1")
    c_local __ctools_plugin_hdfe_markout = ("`ctools_plugin_hdfe_markout'" == "1")
    c_local __ctools_plugin_hdfe_fvlevels = ("`ctools_plugin_hdfe_fvlevels'" == "1")
    capture confirm number 0
end
* ctools-helpers:end
