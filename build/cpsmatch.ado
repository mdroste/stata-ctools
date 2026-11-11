*! version 1.2.6 20261111 github.com/mdroste/stata-ctools

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
         ATE AI(integer 0) POPulation ALTVariance WARNings ///
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
    tempvar order
    quietly gen long `order' = _n

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

    local inputs `varlist' `pscore' `outcome' `mahalanobis'
    foreach input of local inputs {
        if inlist("`input'", "_weight", "_id", "_support", "_pscore", "_treated", "_nn", "_pdif") {
            di as error "input variable `input' conflicts with a matching output"
            exit 111
        }
    }

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
        di as error "Option -ai- only allowed when doing nearest neighbor matching on the covariates (Mahalanobis)"
        exit 198
    }
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

    * psmatch2's generated variables
    foreach v in _treated _support _weight _pscore _id _nn _pdif {
        capture drop `v'
    }
    forvalues n = 1/`neighbor' {
        capture drop _n`n'
    }
    foreach v of local outcome {
        capture drop _`v'
    }

    * Matching sample as psmatch2 builds it: treatment, covariates, the
    * outcomes and the Mahalanobis variables nonmissing, restricted to the
    * score model's e(sample), then a nonmissing score.
    capture markout `touse' `outcome'
    if "`mahalanobis'" != "" {
        fvrevar `mahalanobis', list
        markout `touse' `r(varlist)'
    }

    local score
    if `k_covariates' > 0 {
        local model = cond("`logit'" != "", "logit", "probit")
        `quietly' `model' `depvar' `treatvar' if `touse', nolog
        quietly replace `touse' = e(sample)
        tempvar score
        quietly predict double `score', `index'
        quietly gen double _pscore = `score' if `touse'
        label variable _pscore "cpsmatch: Propensity Score"
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

    * As psmatch2: the treated weigh 1 (0 with ate, then their use as
    * matches), the matched outcome starts at 0 for the rows to be matched.
    quietly gen double _weight = _treated if _support == 1
    char _weight[Type] "aweight"
    if "`ate'" != "" quietly replace _weight = 0 if _treated == 1 & _support == 1
    label variable _weight "cpsmatch: weight of matched controls"
    foreach v of local outcome {
        if "`ate'" == "" quietly gen double _`v' = 0 if _support == 1 & _treated == 1
        else quietly gen double _`v' = 0 if _support == 1
        label variable _`v' "cpsmatch: value of `v' of match(es)"
    }

    if "`warnings'" != "" & "`metric'" == "pscore" {
        _cpsmatch_warn, method(`method') `ate'
    }

    * The plugin matches on the propensity score; Mata matches on the
    * Mahalanobis metric and does local linear regression matching.
    local engine = cond("`metric'" == "pscore" & inlist("`method'", "neighbor", "kernel"), "plugin", "mata")
    local nearest = ("`method'" == "neighbor")

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
    if `nearest' | "`engine'" == "plugin" {
        quietly gen long _nn = 0 if _support == 1
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
        tempvar presup
        quietly gen byte `presup' = _support

        * Several outcomes share one match (it depends on the score alone):
        * the plugin matches once per outcome and stores its matched values.
        local outcomes `outcome'
        if "`outcomes'" == "" local outcomes _cpsmatch_no_outcome_
        local k_outcome = 0
        foreach cur_outcome of local outcomes {
            if "`cur_outcome'" == "_cpsmatch_no_outcome_" local cur_outcome
            local ++k_outcome
            local cur_match : word `k_outcome' of `matchvars'
            * Plugin variables: treatment _pscore presupport _support _weight
            * _nn [_id _n1.._nK _pdif] [outcome _<outcome> [_self_<outcome>]]
            * _treated
            local allvars "`depvar' _pscore `presup' _support _weight _nn"
            local flags "presup=3 nnvar=6"
            local id_idx = 0
            local p = 6
            if `nearest' {
                local id_idx = 7
                local allvars "`allvars' _id `n1'-`nk' _pdif"
                local flags "`flags' ids=8 pdif=`=8 + `neighbor''"
                local p = 8 + `neighbor'
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
            * psmatch2 squares the weights in the default storage type
            if `plugin_stats' & "`c(type)'" == "float" local flags "`flags' w2float"
            if "`ate'" != "" local flags "`flags' ate"
            * The arms are the stored _treated (a code it cannot hold is
            * missing there); the treatment itself marks the sample rows.
            local allvars "`allvars' _treated"
            local flags "`flags' arm=`=`p' + 1'"
            * Format: treat_idx pscore_idx outcome_idx weight_idx id_idx
            *         support_idx nobs method neighbor caliper common kernel
            *         bwidth with_replace ties desc, then flags. The sample
            *         marker goes last, in place of if `touse'.
            local args "1 2 `outcome_idx' 5 `id_idx' 4"
            local args "`args' `nobs' `method_code' `neighbor' `caliper' 0"
            local args "`args' `kernelcode' `bwidth' `with_replace' `handle_ties' `desc_order'"
            plugin call ctools_plugin `allvars' `touse', "cpsmatch `threads_code' `args' `flags' touse"
        }
        * Rows outside the sample follow it in psmatch2's order
        if `nearest' {
            quietly replace _id = `nobs' + sum(!`touse') if !`touse'
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
        forvalues n = 1/`neighbor' {
            quietly recast `idtype' _n`n'
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
    foreach cur_outcome of local outcome {
        local ++k_outcome
        local att = .
        local att_se = .
        local att_t = .
        local atu = .
        local ate_v = .
        local seatu = .
        local seate = .
        if `plugin_stats' {
            * The plugin's statistics of the last outcome's call; rerun for
            * earlier outcomes, whose scalars it has replaced.
            if `k_outcome' < `: word count `outcome'' {
                _cpsmatch_effects `cur_outcome', neighbor(`neighbor')
                local y1_bar = r(m1t)
                local y0_bar = r(m0t)
                local att = r(att)
                local att_se = r(seatt)
                local y1_unmatched = r(u1u)
                local y0_unmatched = r(u0u)
                local unmatched_se = r(seols)
            }
            else {
                local y1_bar = scalar(_cpsmatch_y1)
                local n1 = scalar(_cpsmatch_n1)
                local var1 = scalar(_cpsmatch_var1)
                local y0_bar = scalar(_cpsmatch_y0)
                local n0 = scalar(_cpsmatch_n0)
                if `n1' > 0 & `n0' > 0 {
                    local att = `y1_bar' - `y0_bar'
                    * psmatch2 default ATT variance: unweighted group variances,
                    * with squared match-frequency weights for reused controls.
                    local var0 = scalar(_cpsmatch_var0)
                    local sum_weight_squared = scalar(_cpsmatch_w2)
                    local att_se = sqrt(`var1'/`n1' + `var0'*`sum_weight_squared'/`n1'^2)
                }
                local y1_unmatched = scalar(_cpsmatch_y1_u)
                local y0_unmatched = scalar(_cpsmatch_y0_u)
                * psmatch2 displays the OLS SE from regressing outcome on
                * treatment: its pooled residual variance, where a
                * singleton group contributes zero SSE.
                local var1_u = scalar(_cpsmatch_var1_u)
                local n1_u = scalar(_cpsmatch_n1_u)
                local var0_u = scalar(_cpsmatch_var0_u)
                local n0_u = scalar(_cpsmatch_n0_u)
                local sse1_u = cond(`n1_u' > 1, (`n1_u'-1)*`var1_u', 0)
                local sse0_u = cond(`n0_u' > 1, (`n0_u'-1)*`var0_u', 0)
                local pooled_var = (`sse1_u' + `sse0_u')/(`n1_u' + `n0_u' - 2)
                local unmatched_se = sqrt(`pooled_var'*(1/`n1_u' + 1/`n0_u'))
            }
        }
        else {
            local skipse = ""
            if `smooth' == 2 | ("`method'" == "llr" & "`kerneltype'" == "tricube") local skipse skipse
            _cpsmatch_effects `cur_outcome', neighbor(`neighbor') ai(`ai') `ate' ///
                `population' `altvariance' `skipse'
            local y1_bar = r(m1t)
            local y0_bar = r(m0t)
            local att = r(att)
            local att_se = r(seatt)
            local y1_unmatched = r(u1u)
            local y0_unmatched = r(u0u)
            local unmatched_se = r(seols)
            if "`ate'" != "" {
                local m1u = r(m1u)
                local m0u = r(m0u)
                local atu = r(atu)
                local ate_v = r(ate)
                local seatu = r(seatu)
                local seate = r(seate)
            }
        }
        if !missing(`att_se') & `att_se' > 0 local att_t = `att' / `att_se'

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
        local last_att_t = `att_t'
    }
    if `shown_header' {
        if `ai' == 0 {
            if `last_seatt' != . di as text "Note: S.E. does not take into account that the propensity score is estimated."
        }
        else if "`population'" == "" di as text "Note: Sample S.E."
        else di as text "Note: Population S.E."
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
    syntax varname, [neighbor(integer 1) ai(integer 0) ate population altvariance skipse]
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
            * with altvariance. As psmatch2, the terms are generated in the
            * default storage type.
            tempvar w shat vt vu ve
            generate `w' = max(_weight, 0)
            if "`altvariance'" == "" generate `shat' = (`ai'/(`ai' + 1))*(`v' - _self_`v')^2
            else generate `shat' = _self_`v'
            if "`population'" == "" {
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
            tempvar w2
            summarize `v' if _treated == 1 & _support == 1
            local var1 = r(Var)
            summarize `v' if _treated == 0 & _weight < .
            local var0 = r(Var)
            generate `w2' = _weight^2 if _treated == 0
            summarize `w2', meanonly
            local wtot = r(sum)
            if "`skipse'" != "" scalar `seatt' = .
            else scalar `seatt' = sqrt(`var1'/`n1' + `var0'*`wtot'/`n1'^2)
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
        quietly _rmdcoll _treated `mahalanobis' if `touse', expand
        local mahalanobis `r(varlist)'
        fvrevar `mahalanobis'
        local xvars `r(varlist)'
        if "`w'" == "" {
            tempname XX0 XX1
            quietly matrix accum `XX0' = `xvars' if _treated == 0, deviations noconstant
            quietly matrix accum `XX1' = `xvars' if _treated == 1, deviations noconstant
            quietly count if _treated <= 1 & _support == 1
            matrix `W' = syminv((`XX0' + `XX1')/(r(N) - 2))
        }
        else {
            matrix `W' = `w'
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
    real matrix X, W, IDS, Y, M, d, best
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
    q = J(rows(X), 1, .)
    for (i = first; i <= last; i++) q[i] = X[i, .] * W * X[i, .]'

    for (i = t0; i <= t1; i++) {
        if (self & !ate & weight[i] == 0) continue
        d = q[|c0 \ c1|] - 2 * X[|c0, 1 \ c1, .|] * (W * X[i, .]') :+ q[i]
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
    real matrix X, W, Y, M
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
    if (mahal) {
        q = J(rows(X), 1, .)
        for (i = 1; i <= rows(X); i++) q[i] = X[i, .] * W * X[i, .]'
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
        if (mahal) d = q[donors] - 2 * X[donors, .] * (W * X[i, .]') :+ q[i]
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

end
