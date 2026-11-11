*! version 1.2.6 20261111 github.com/mdroste/stata-ctools

program define crangestat, rclass
    version 14.1
    * A failed or interrupted call leaves the data as they were without a
    * snapshot of the dataset: crangestat only adds its result variables at
    * the end, and the plugin writes into nothing else, so _ctools_undo drops
    * them and resets c(changed). Option snapshot preserves the whole dataset
    * instead (on disk above c(max_preservemem)); option fast undoes nothing.
    _ctools_mode `0'
    * Break is recognized only inside _crangestat_impl.
    nobreak {
        if "`mode'" == "snapshot" preserve
        else if "`mode'" == "" _ctools_undo record CTOOLS_u_crangestat_
        capture noisily break _crangestat_impl `0'
        local rc = _rc
        if `rc' {
            if "`mode'" == "snapshot" restore
            else if "`mode'" == "" _ctools_undo rollback CTOOLS_u_crangestat_
            exit `rc'
        }
        if "`mode'" == "snapshot" restore, not
        else if "`mode'" == "" _ctools_undo clear CTOOLS_u_crangestat_
        * rangestat's results: the new variables, in the order of the
        * statistics, also in the caller's local(), and rangestat's version
        if "`__crs_local'" != "" c_local `__crs_local' `__crs_newvars'
        return local newvars `__crs_newvars'
        return local rangestat_version 1.1.1
    }
end

program define _crangestat_impl
    version 14.1

    * Check observation limit (Stata plugin API limitation)
    if _N > 2147483647 {
        di as error "ctools does not support datasets exceeding 2^31 (2.147 billion) observations"
        di as error "This is a limitation of Stata's plugin API"
        exit 920
    }

    * Parse the syntax as rangestat does:
    * crangestat (stat) [newvar=]varname [...] [if] [in], interval(keyvar low high) [by(varlist) excludeself verbose threads(#)]
    syntax anything(name=stat_part id="statistics" equalok) [if] [in], ///
        Interval(string) [BY(varlist) EXCLudeself CASEwise Describe local(name) noCHeck ///
        Verbose THReads(integer 0) SNAPshot FAST]

    * Freeze incoming-r() qualifiers before parsing helpers replace r().
    if `"`if'"' != "" {
        tempvar __incoming_sample
        mark `__incoming_sample' `if' `in'
        local if "if `__incoming_sample'"
    }
    if !inrange(`threads', 0, 256) {
        di as error "threads() must be between 0 and 256"
        exit 198
    }
    * nocheck skips rangestat's shortcut for observations with the same
    * interval bounds, which only saves time there: it changes no result.

    * Parse interval specification: keyvar low high, separated by spaces or
    * commas as in rangestat, which refuses anything more
    local interval : subinstr local interval "," " ", all
    tokenize "`interval'"
    local keyvar "`1'"
    local int_low "`2'"
    local int_high "`3'"
    if `"`4'"' != "" {
        di as error "extra argument in interval() option: `4'"
        exit 198
    }

    if "`keyvar'" == "" | "`int_low'" == "" | "`int_high'" == "" {
        di as error "crangestat: interval() requires keyvar low high"
        exit 198
    }

    * Validate key variable exists and is numeric
    capture confirm variable `keyvar'
    if _rc != 0 {
        di as error "crangestat: variable `keyvar' not found"
        exit 111
    }
    capture confirm numeric variable `keyvar'
    if _rc != 0 {
        di as error "crangestat: key variable `keyvar' must be numeric"
        exit 109
    }

    * Convert interval bounds (handle . for +/- infinity). As in rangestat, a
    * bound may be a numeric variable: the bound itself (not an offset),
    * with a missing value unbounded. Then both bounds go to the plugin as
    * variables, per observation.
    local low_is_missing = 0
    local high_is_missing = 0
    local var_bounds = 0
    foreach __b in low high {
        if regexm("`int_`__b''", "^[.][a-z]?$") local int_`__b' .
        capture confirm numeric variable `int_`__b''
        if !_rc {
            unab int_`__b' : `int_`__b'', max(1)
            local var_bounds = 1
        }
        else {
            capture confirm number `int_`__b''
            if !_rc {
                if missing(`int_`__b'') local int_`__b' .
            }
        }
    }
    if `var_bounds' {
        tempvar bound_low bound_high
        foreach __b in low high {
            capture confirm numeric variable `int_`__b''
            if !_rc {
                quietly generate double `bound_`__b'' = cond(`int_`__b'' > ., ., `int_`__b'')
            }
            else if "`int_`__b''" == "." {
                quietly generate double `bound_`__b'' = .
            }
            else {
                capture confirm number `int_`__b''
                if _rc {
                    di as error "was expecting a numeric variable, a number, or a system missing value for the interval `__b': `int_`__b''"
                    exit 198
                }
                quietly generate double `bound_`__b'' = `keyvar' + `int_`__b''
            }
        }
        local int_low "0"
        local int_high "0"
    }

    if "`int_low'" == "." {
        local low_is_missing = 1
        local int_low_val = 0
    }
    else {
        capture confirm number `int_low'
        if _rc != 0 {
            di as error "crangestat: invalid interval low bound '`int_low''"
            exit 198
        }
        local int_low_val = `int_low'
    }

    if "`int_high'" == "." {
        local high_is_missing = 1
        local int_high_val = 0
    }
    else {
        capture confirm number `int_high'
        if _rc != 0 {
            di as error "crangestat: invalid interval high bound '`int_high''"
            exit 198
        }
        local int_high_val = `int_high'
    }

    * Start timing
    local __do_timing = ("`verbose'" != "")
    if `__do_timing' {
        local __ct_elapsed90 = 0
        local __ct_elapsed91 = 0
        local __ct_elapsed92 = 0
        if "`__ct_elapsed90'" == "" local __ct_elapsed90 = 0
        quietly _ctools_clock __ct_start90
        if "`__ct_elapsed91'" == "" local __ct_elapsed91 = 0
        quietly _ctools_clock __ct_start91
    }

    * Parse stat specifications
    * Format: (stat) [newvar=]varname [(stat) [newvar=]varname ...]
    local stat_codes ""
    local source_vars ""
    local result_vars ""
    local source_var_map ""
    local nstats = 0
    local nsource = 0
    local source_var_list ""

    local remaining "`stat_part'"
    local nflex = 0
    local flex_results
    local ncw = 0
    local cw_fix
    local cw_hidden
    local newvars_order
    * As in rangestat, a specification without a statistic means (mean).
    if substr(strtrim("`remaining'"), 1, 1) != "(" local remaining "(mean) `remaining'"

    while "`remaining'" != "" {
        * Skip leading whitespace
        local remaining = strtrim("`remaining'")
        if "`remaining'" == "" continue, break

        * Expect opening parenthesis
        if substr("`remaining'", 1, 1) != "(" {
            di as error "crangestat: expected '(' before statistic name"
            exit 198
        }

        * Find closing parenthesis
        local paren_end = strpos("`remaining'", ")")
        if `paren_end' == 0 {
            di as error "crangestat: missing ')' after statistic name"
            exit 198
        }

        * Extract stat name
        local stat_name = strtrim(substr("`remaining'", 2, `paren_end' - 2))
        local remaining = strtrim(substr("`remaining'", `paren_end' + 1, .))

        * Validate stat name and get code
        local stat_code = -1
        if "`stat_name'" == "count" local stat_code = 0
        else if "`stat_name'" == "mean" local stat_code = 1
        else if "`stat_name'" == "sum" local stat_code = 2
        else if "`stat_name'" == "min" local stat_code = 3
        else if "`stat_name'" == "max" local stat_code = 4
        else if "`stat_name'" == "sd" local stat_code = 5
        else if "`stat_name'" == "variance" local stat_code = 6
        else if "`stat_name'" == "median" local stat_code = 7
        else if "`stat_name'" == "iqr" local stat_code = 8
        else if "`stat_name'" == "first" local stat_code = 9
        else if "`stat_name'" == "last" local stat_code = 10
        else if "`stat_name'" == "firstnm" local stat_code = 11
        else if "`stat_name'" == "lastnm" local stat_code = 12
        else if "`stat_name'" == "p1" local stat_code = 13
        else if "`stat_name'" == "p5" local stat_code = 14
        else if "`stat_name'" == "p10" local stat_code = 15
        else if "`stat_name'" == "p25" local stat_code = 16
        else if "`stat_name'" == "p75" local stat_code = 17
        else if "`stat_name'" == "p90" local stat_code = 18
        else if "`stat_name'" == "p95" local stat_code = 19
        else if "`stat_name'" == "p99" local stat_code = 20
        else if "`stat_name'" == "skewness" local stat_code = 21
        else if "`stat_name'" == "kurtosis" local stat_code = 22
        else if "`stat_name'" == "obs" local stat_code = 23
        else if "`stat_name'" == "missing" local stat_code = 24

        * rangestat's flexible statistics: (cov), (corr), (reg) and a
        * user-written Mata function, computed in Mata over each window
        local is_flex = 0
        if `stat_code' == -1 {
            if inlist("`stat_name'", "cov", "corr", "reg") local is_flex = 1
            else if regexm("`stat_name'", "^[A-Za-z_][A-Za-z0-9_]*$") {
                mata: st_local("is_flex", strofreal(findexternal("`stat_name'()") != NULL))
            }
        }
        if `stat_code' == -1 & !`is_flex' {
            di as error "crangestat: unknown statistic '`stat_name''"
            exit 3499
        }
        if `is_flex' {
            local next_paren = strpos("`remaining'", "(")
            if `next_paren' == 0 {
                local var_spec = strtrim("`remaining'")
                local remaining ""
            }
            else {
                local var_spec = strtrim(substr("`remaining'", 1, `next_paren' - 1))
                local remaining = substr("`remaining'", `next_paren', .)
            }
            if "`var_spec'" == "" {
                di as error "crangestat: missing variable after (`stat_name')"
                exit 198
            }
            * A variable range may have spaces around the dash, as in rangestat
            local var_spec = ustrregexra(`"`var_spec'"', "\s*-\s*", "-")
            unab var_spec : `var_spec'
            confirm numeric variable `var_spec'
            local nflexvars : word count `var_spec'
            if inlist("`stat_name'", "cov", "corr") & `nflexvars' != 2 {
                di as error "`stat_name' stat requires exactly 2 variables; you included: `var_spec'"
                exit 198
            }
            if "`stat_name'" == "reg" & `nflexvars' < 2 {
                di as error "reg stat requires a dependent variable and at least one predictor: `var_spec'"
                exit 198
            }
            local ++nflex
            local flex_name_`nflex' `stat_name'
            local flex_vars_`nflex' `var_spec'
            if "`stat_name'" == "cov" local flex_out_`nflex' cov_nobs cov_x
            else if "`stat_name'" == "corr" local flex_out_`nflex' corr_nobs corr_x
            else if "`stat_name'" == "reg" {
                gettoken __dep __indep : var_spec
                local flex_out_`nflex' reg_nobs reg_r2 reg_adj_r2
                foreach __v in `__indep' cons {
                    local flex_out_`nflex' `flex_out_`nflex'' b_`__v'
                }
                foreach __v in `__indep' cons {
                    local flex_out_`nflex' `flex_out_`nflex'' se_`__v'
                }
            }
            local flex_results `flex_results' `flex_out_`nflex''
            local newvars_order `newvars_order' @flex`nflex'
            continue
        }

        * Parse variable specification: [newvar=]varname
        * Find next ( or end
        local next_paren = strpos("`remaining'", "(")
        if `next_paren' == 0 {
            local var_spec = strtrim("`remaining'")
            local remaining ""
        }
        else {
            local var_spec = strtrim(substr("`remaining'", 1, `next_paren' - 1))
            local remaining = substr("`remaining'", `next_paren', .)
        }

        if "`var_spec'" == "" {
            di as error "crangestat: missing variable after (`stat_name')"
            exit 198
        }

        * Each specification holds [newvar=]varname pairs and varlists, as in
        * rangestat; a varlist gives one result per variable, named
        * varname_statname.
        local var_spec = ustrregexra(`"`var_spec'"', "\s*-\s*", "-")
        local var_spec : subinstr local var_spec "=" " = ", all
        local spec_sources ""
        local spec_results ""
        while `"`var_spec'"' != "" {
            gettoken word var_spec : var_spec
            gettoken eq : var_spec
            if `"`word'"' == "=" {
                di as error "crangestat: missing new variable name in (`stat_name')"
                exit 102
            }
            if `"`eq'"' == "=" {
                gettoken eq var_spec : var_spec
                gettoken source_var var_spec : var_spec
                if "`source_var'" == "" | "`source_var'" == "=" {
                    di as error "crangestat: missing variable after `word'= in (`stat_name')"
                    exit 198
                }
                capture unab source_var : `source_var', max(1)
                if _rc {
                    di as error "crangestat: variable `source_var' not found"
                    exit 111
                }
                local spec_sources `spec_sources' `source_var'
                local spec_results `spec_results' `word'
            }
            else {
                capture unab expanded : `word'
                if _rc {
                    di as error "crangestat: variable `word' not found"
                    exit 111
                }
                foreach source_var of local expanded {
                    local spec_sources `spec_sources' `source_var'
                    local spec_results `spec_results' `source_var'_`stat_name'
                }
            }
        }

        local nspec : word count `spec_sources'
        * casewise: within each window, rangestat drops the observations
        * with a missing value of any variable of this statistic. The plugin
        * gets copies of them that are missing on those observations, and a
        * count of the complete observations of each window (a hidden
        * statistic) makes the results of a window without any missing. On the
        * complete observations, first and last are firstnm and lastnm, obs
        * is that count and missing is 0.
        local spec_labels `spec_sources'
        local cw_n
        if "`casewise'" != "" {
            foreach source_var of local spec_sources {
                capture confirm numeric variable `source_var'
                if _rc {
                    di as error "crangestat: source variable `source_var' must be numeric"
                    exit 109
                }
            }
            local ++ncw
            tempvar cw_complete cw_n
            quietly generate byte `cw_complete' = 1
            foreach source_var of local spec_sources {
                quietly replace `cw_complete' = . if missing(`source_var')
            }
            local cw_masked
            foreach source_var of local spec_sources {
                tempvar cw_copy
                quietly generate double `cw_copy' = `source_var' if `cw_complete' == 1
                local cw_masked `cw_masked' `cw_copy'
            }
            local spec_sources `cw_masked'
        }
        forvalues spec_j = 1/`nspec' {
        local source_var : word `spec_j' of `spec_sources'
        local result_var : word `spec_j' of `spec_results'
        local source_label : word `spec_j' of `spec_labels'
        local stat_code_j = `stat_code'
        if "`cw_n'" != "" {
            if inlist(`stat_code', 23, 24) {
                local source_var `cw_complete'
                local stat_code_j = 0
            }
            else if `stat_code' == 9 local stat_code_j = 11
            else if `stat_code' == 10 local stat_code_j = 12
            local cw_fix `cw_fix' `result_var' `stat_code' `cw_n'
        }

        capture confirm numeric variable `source_var'
        if _rc != 0 {
            di as error "crangestat: source variable `source_var' must be numeric"
            exit 109
        }
        confirm name `result_var'

        * Destination creation is deferred until the entire request is valid.

        * Track source variable index
        local src_idx = 0
        local found = 0
        foreach sv of local source_var_list {
            local src_idx = `src_idx' + 1
            if "`sv'" == "`source_var'" {
                local found = 1
                continue, break
            }
        }
        if `found' == 0 {
            local source_var_list "`source_var_list' `source_var'"
            local nsource = `nsource' + 1
            local src_idx = `nsource'
        }

        * Store mapping (0-based for C)
        local source_var_map "`source_var_map' `=`src_idx' - 1'"

        local stat_codes "`stat_codes' `stat_code_j'"
        local result_vars "`result_vars' `result_var'"
        local result_label_`=`nstats'+1' "`stat_name' of `source_label'"
        local nstats = `nstats' + 1
        local newvars_order `newvars_order' `result_var'
        }
        * The hidden count of complete observations for casewise
        if "`cw_n'" != "" {
            local src_idx = 0
            local found = 0
            foreach sv of local source_var_list {
                local ++src_idx
                if "`sv'" == "`cw_complete'" {
                    local found = 1
                    continue, break
                }
            }
            if !`found' {
                local source_var_list "`source_var_list' `cw_complete'"
                local ++nsource
                local src_idx = `nsource'
            }
            local source_var_map "`source_var_map' `=`src_idx' - 1'"
            local stat_codes "`stat_codes' 0"
            local result_vars "`result_vars' `cw_n'"
            local result_label_`=`nstats'+1' "complete observations"
            local nstats = `nstats' + 1
            local cw_hidden `cw_hidden' `cw_n'
        }
    }

    if `nstats' == 0 & `nflex' == 0 {
        di as error "crangestat: no statistics specified"
        exit 198
    }
    * rangestat rejects a result name used twice, or already defined, with r(110).
    local seen_results ""
    foreach output in `result_vars' `flex_results' {
        if `: list output in seen_results' {
            di as error "cannot create -`output'-; variable(s) already defined"
            exit 110
        }
        local seen_results `seen_results' `output'
    }
    local visible_results : list result_vars - cw_hidden
    if "`visible_results'`flex_results'" != "" _ctools_newvars `visible_results' `flex_results'

    * The sample follows rangestat: if/in and a nonmissing key. Excluded rows
    * take a missing key, so they neither join a window nor get a result.
    if `"`if'`in'"' != "" {
        marksample touse, novarlist
        tempvar sample_key
        quietly generate double `sample_key' = `keyvar' if `touse'
        local keyvar `sample_key'
    }
    quietly count if !missing(`keyvar')
    if r(N) == 0 error 2000
    local result_index=0
    foreach output of local result_vars {
        local ++result_index
        quietly generate double `output' = .
        label variable `output' "`result_label_`result_index''"
    }
    * Validation hook: fail as Break would once the result variables exist
    if "$CTOOLS_crangestat_fault" == "outputs" error 1

    * Load plugin
    _ctools_load
    * Stata scopes plugin registrations to the calling ado program.
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    capture confirm number 0

    * strL by() values the plugin cannot read (binary, or over 2,045 bytes)
    * are grouped through stand-in ranks. The Mata statistics below group by
    * the variables themselves.
    local plugin_by `by'
    if "`by'" != "" & `nstats' > 0 {
        local stubs
        foreach v of local by {
            tempvar stub
            local stubs `stubs' `stub'
        }
        _ctools_strlkey `by', generate(`stubs')
        local plugin_by `__ctools_keys'
    }

    * Get variable indices
    * Build list of ALL dataset variables - this determines plugin varlist order
    * The plugin interface uses 1-based indices relative to the varlist passed to it
    unab all_vars : *

    * Indices are dataset positions (the plugin receives every variable in
    * order). Look up only the variables used: a local named after each
    * dataset variable is an invalid name once a name has 24+ characters.
    _crangestat_varpos key_idx : `keyvar'
    * Variable bounds: two more sources after the statistics' ones
    if `var_bounds' {
        local source_var_list `source_var_list' `bound_low' `bound_high'
        local nsource = `nsource' + 2
    }
    _crangestat_varpos source_indices : `source_var_list'
    _crangestat_varpos result_indices : `result_vars'

    * By-variable indices
    local nby = 0
    local by_indices ""
    if "`by'" != "" {
        local nby : word count `plugin_by'
        _crangestat_varpos by_indices : `plugin_by'
    }

    * Build options string
    local opts ""
    if "`excludeself'" != "" local opts "`opts' excludeself"
    if "`verbose'" != "" local opts "`opts' verbose"
    local opts "`opts' low=`int_low_val' high=`int_high_val'"
    if `low_is_missing' local opts "`opts' low_miss"
    if `high_is_missing' local opts "`opts' high_miss"
    if `var_bounds' local opts "`opts' varbounds lowsrc=`=`nsource'-2' highsrc=`=`nsource'-1'"

    local threads_opt ""
    if `threads' > 0 {
        local threads_opt "threads(`threads')"
    }

    * Build argument string for plugin
    * Format: nstats nsource nby key_idx source_indices... result_indices... by_indices... stat_types... source_var_map... [options]
    local plugin_args "`nstats' `nsource' `nby' `key_idx' `source_indices' `result_indices' `by_indices' `stat_codes' `source_var_map' `opts'"

    * Call plugin
    if `__do_timing' {
        quietly _ctools_clock __ct_now
        local __ct_elapsed91 = `__ct_elapsed91' + (`__ct_now' - `__ct_start91')
        if "`__ct_elapsed92'" == "" local __ct_elapsed92 = 0
        quietly _ctools_clock __ct_start92
    }
    if `nstats' > 0 plugin call ctools_plugin `all_vars', "crangestat `threads_opt' `plugin_args'"
    * casewise: a window without complete observations has no result (obs
    * is their count; missing is 0 when there are some)
    while "`cw_fix'" != "" {
        gettoken cw_r cw_fix : cw_fix
        gettoken cw_code cw_fix : cw_fix
        gettoken cw_count cw_fix : cw_fix
        quietly replace `cw_r' = . if `cw_count' == 0 | missing(`cw_count')
        if `cw_code' == 24 quietly replace `cw_r' = 0 if !missing(`cw_r')
    }

    * Flexible statistics, in Mata over the same windows: the bounds as
    * rangestat forms them (a number is added to the key; a variable is the
    * bound; a missing bound is unbounded) and its by() groups
    if `nflex' > 0 {
        tempvar flex_low flex_high flex_group flex_touse
        foreach __b in low high {
            if `var_bounds' quietly generate double `flex_`__b'' = `bound_`__b''
            else if ``__b'_is_missing' quietly generate double `flex_`__b'' = .
            else quietly generate double `flex_`__b'' = `keyvar' + `int_`__b'_val'
        }
        if "`by'" != "" quietly egen long `flex_group' = group(`by'), missing
        else quietly generate byte `flex_group' = 1
        quietly generate byte `flex_touse' = !missing(`keyvar')
        forvalues f = 1/`nflex' {
            local fname `flex_name_`f''
            local fvars `flex_vars_`f''
            local builtin = inlist("`fname'", "cov", "corr", "reg")
            local cw = `builtin' | "`casewise'" != ""
            local xs = "`excludeself'" != ""
            mata: _crangestat_flex("`keyvar'", "`flex_low'", "`flex_high'", "`flex_group'", ///
                "`flex_touse'", "`fvars'", "`fname'", `builtin', `cw', `xs', "flex_tmp")
            local ntmp : word count `flex_tmp'
            if `builtin' {
                local want : word count `flex_out_`f''
                if `ntmp' != `want' {
                    di as error "no result for all obs: `fname' `fvars'"
                    drop `flex_tmp'
                    continue
                }
                local flexnames_`f' `flex_out_`f''
                local reg_label "regress `fvars'"
                if strlen("`reg_label'") > 50 local reg_label = substr("`reg_label'", 1, 50) + "..."
                local nb = (`ntmp' - 3) / 2
                forvalues j = 1/`ntmp' {
                    local from : word `j' of `flex_tmp'
                    local to : word `j' of `flex_out_`f''
                    rename `from' `to'
                    if "`fname'" == "cov" {
                        if `j' == 1 label variable `to' "number of observations - cov"
                        else label variable `to' "covariance of `fvars'"
                    }
                    else if "`fname'" == "corr" {
                        if `j' == 1 label variable `to' "number of observations - corr `fvars'"
                        else label variable `to' "correlation of `fvars'"
                    }
                    else {
                        if `j' == 1 label variable `to' "number of obs - regress `reg_label'"
                        else if `j' == 2 label variable `to' "R-squared - `reg_label'"
                        else if `j' == 3 label variable `to' "adj. R-squared - `reg_label'"
                        else if `j' <= 3 + `nb' label variable `to' "coef - `reg_label'"
                        else label variable `to' "std err - `reg_label'"
                    }
                }
            }
            else {
                forvalues j = 1/`ntmp' {
                    local from : word `j' of `flex_tmp'
                    confirm new variable `fname'`j'
                    rename `from' `fname'`j'
                    label variable `fname'`j' "`fname' of `fvars'"
                    local flexnames_`f' `flexnames_`f'' `fname'`j'
                }
            }
        }
    }
    * The new variables in the order of the statistics, as rangestat
    * returns them; describe and local() as there
    local newvars
    foreach w of local newvars_order {
        if substr("`w'", 1, 5) == "@flex" {
            local f = substr("`w'", 6, .)
            local newvars `newvars' `flexnames_`f''
        }
        else local newvars `newvars' `w'
    }
    * rangestat adds them in that order
    if "`newvars'" != "" order `newvars', last
    if "`describe'" != "" & "`newvars'" != "" describe `newvars'
    c_local __crs_newvars `newvars'
    c_local __crs_local `local'
    * Validation hook: fail after the plugin filled the result variables
    if "$CTOOLS_crangestat_fault" == "post" error 1
    if `__do_timing' {
        quietly _ctools_clock __ct_now
        local __ct_elapsed92 = `__ct_elapsed92' + (`__ct_now' - `__ct_start92')
    }

    * Display summary or timing
    if `__do_timing' {
        quietly _ctools_clock __ct_now
        local __ct_elapsed90 = `__ct_elapsed90' + (`__ct_now' - `__ct_start90')
        * Elapsed seconds are held in caller-local timing state.
        local __time_total = `__ct_elapsed90'

        * Extract sub-timer values
        * Elapsed seconds are held in caller-local timing state.
        local __time_preplugin = `__ct_elapsed91'
        * Elapsed seconds are held in caller-local timing state.
        local __time_plugin = `__ct_elapsed92'

        * Calculate plugin call overhead
        local __plugin_time_total = 0
        if `nstats' > 0 local __plugin_time_total = _crangestat_time_total
        local __plugin_call_overhead = `__time_plugin' - `__plugin_time_total'

        di as text ""
        di as text "{hline 55}"
        di as text "crangestat timing breakdown:"
        di as text "{hline 55}"
        if `nstats' > 0 {
        di as text "  C plugin internals:"
        di as text "    Data load:              " as result %8.4f _crangestat_time_load " sec"
        di as text "    Sort:                   " as result %8.4f _crangestat_time_sort " sec"
        di as text "    Group detection:        " as result %8.4f _crangestat_time_groups " sec"
        di as text "    Compute statistics:     " as result %8.4f _crangestat_time_compute " sec"
        di as text "    Store results:          " as result %8.4f _crangestat_time_store " sec"
        di as text "  {hline 53}"
        di as text "    C plugin total:         " as result %8.4f _crangestat_time_total " sec"
        di as text ""
        }
        di as text "  Stata overhead:"
        di as text "    Pre-plugin setup:       " as result %8.4f `__time_preplugin' " sec"
        di as text "    Plugin call overhead:   " as result %8.4f `__plugin_call_overhead' " sec"
        di as text "  {hline 53}"
        local __stata_overhead = `__time_preplugin' + `__plugin_call_overhead'
        di as text "    Stata overhead total:   " as result %8.4f `__stata_overhead' " sec"
        di as text "{hline 55}"
        di as text "    Wall clock total:       " as result %8.4f `__time_total' " sec"
        di as text "{hline 55}"

        * Display thread diagnostics
        local __threads_max = .
        if `nstats' > 0 local __threads_max = _crangestat_threads_max
        if !missing(`__threads_max') {
            capture local __openmp_enabled = _crangestat_openmp_enabled
            if _rc != 0 local __openmp_enabled = 0
            di as text ""
            di as text "  Thread diagnostics:"
            di as text "    OpenMP enabled:         " as result %8.0f `__openmp_enabled'
            di as text "    Max threads available:  " as result %8.0f `__threads_max'
            di as text "{hline 55}"
        }

        * Clean up timing scalars
        capture scalar drop _crangestat_time_load _crangestat_time_sort _crangestat_time_groups
        capture scalar drop _crangestat_time_compute _crangestat_time_store _crangestat_time_total
        capture scalar drop _crangestat_threads_max _crangestat_openmp_enabled
    }
    else {
        di as text "(`=`nstats' - `ncw' + `nflex'' range statistics computed)"
    }
end

* Dataset positions of the named variables (abbreviations allowed), as a list.
program define _crangestat_varpos
    gettoken macname 0 : 0
    gettoken colon 0 : 0
    local pos ""
    if `"`0'"' != "" {
        mata: st_local("pos", invtokens(strofreal(st_varindex(tokens(st_local("0")), 1))))
    }
    c_local `macname' "`pos'"
end

* The flexible statistics of rangestat, for each observation of the sample:
* the rows of xvars whose key lies in [low, high] within its group (a
* missing low is unbounded below, a missing high above), with the current
* row's values missing under excludeself and rows with a missing value
* dropped under casewise deletion. A nonempty window gives the function's
* row vector; the result columns grow to the widest one. The tempvars that
* hold the results are returned in the local named by outlocal.
capture mata: mata drop _crangestat_flex()
capture mata: mata drop _crangestat_cov()
capture mata: mata drop _crangestat_corr()
capture mata: mata drop _crangestat_reg()
capture mata: mata drop _crangestat_first()
capture mata: mata drop _crangestat_last()
mata:
real rowvector _crangestat_cov(real matrix X)
{
    real matrix V
    V = variance(X)
    return((rows(X), V[1, 2]))
}

real rowvector _crangestat_corr(real matrix X)
{
    real matrix C
    C = correlation(X)
    return((rows(X), C[1, 2]))
}

/* OLS of the first column on the others and a constant: the number of
 * rows, R-squared, adjusted R-squared, coefficients and standard errors.
 * Only the count is returned when the window has no more rows than
 * parameters or a regressor would be omitted. */
real rowvector _crangestat_reg(real matrix Xall)
{
    real colvector y, b, e, se
    real matrix X, XXi
    real scalar n, k, ybar, tss, mss, r2

    y = Xall[., 1]
    n = rows(Xall)
    X = (Xall[., 2..cols(Xall)], J(n, 1, 1))
    k = cols(X)
    if (n <= k) return(n)
    XXi = invsym(quadcross(X, X))
    if (diag0cnt(XXi)) return(n)
    b = XXi * quadcross(X, y)
    e = y - X * b
    se = sqrt(diagonal((quadsum(e :^ 2) / (n - k)) * XXi))
    ybar = mean(y)
    tss = sum((y :- ybar) :^ 2)
    mss = sum((X * b :- ybar) :^ 2)
    r2 = mss / tss
    return((n, r2, 1 - (1 - r2) * (n - 1) / (n - k), b', se'))
}

/* First index in sorted keys K[a..b] with K >= v, or b + 1 */
real scalar _crangestat_first(real colvector K, real scalar a, real scalar b, real scalar v)
{
    real scalar lo, hi, mid
    lo = a
    hi = b + 1
    while (lo < hi) {
        mid = floor((lo + hi) / 2)
        if (K[mid] < v) lo = mid + 1
        else hi = mid
    }
    return(lo)
}

/* Last index in sorted keys K[a..b] with K <= v, or a - 1 */
real scalar _crangestat_last(real colvector K, real scalar a, real scalar b, real scalar v)
{
    real scalar lo, hi, mid
    lo = a
    hi = b + 1
    while (lo < hi) {
        mid = floor((lo + hi) / 2)
        if (K[mid] <= v) lo = mid + 1
        else hi = mid
    }
    return(lo - 1)
}

void _crangestat_flex(string scalar keyvar, string scalar lowvar, string scalar highvar,
    string scalar groupvar, string scalar touse, string scalar xvars, string scalar fname,
    real scalar builtin, real scalar casewise, real scalar exself, string scalar outlocal)
{
    real colvector K, L, H, G, ord, Ks
    real matrix X, R, Xsub, info
    real rowvector res
    real scalar n, g, i, a, b, ws, we, lo, hi, row
    string rowvector names
    pointer(real rowvector function) scalar f

    if (builtin) {
        if (fname == "cov") f = &_crangestat_cov()
        else if (fname == "corr") f = &_crangestat_corr()
        else f = &_crangestat_reg()
    }
    else f = findexternal(fname + "()")

    K = st_data(., keyvar, touse)
    L = st_data(., lowvar, touse)
    H = st_data(., highvar, touse)
    G = st_data(., groupvar, touse)
    X = st_data(., tokens(xvars), touse)
    n = rows(K)
    ord = order((G, K, (1::n)), (1, 2, 3))
    Ks = K[ord]
    info = panelsetup(G[ord], 1)
    R = J(n, 1, .)
    for (g = 1; g <= rows(info); g++) {
        a = info[g, 1]
        b = info[g, 2]
        for (i = a; i <= b; i++) {
            row = ord[i]
            lo = L[row]
            hi = H[row]
            ws = (lo >= .) ? a : _crangestat_first(Ks, a, b, lo)
            we = (hi >= .) ? b : _crangestat_last(Ks, a, b, hi)
            if (we < ws) continue
            Xsub = X[ord[ws..we], .]
            if (exself & i >= ws & i <= we) Xsub[i - ws + 1, .] = J(1, cols(X), .)
            if (casewise) Xsub = select(Xsub, rowmissing(Xsub) :== 0)
            if (!rows(Xsub)) continue
            res = (*f)(Xsub)
            if (cols(res) > cols(R)) R = (R, J(n, cols(res) - cols(R), .))
            if (cols(res)) R[row, 1..cols(res)] = res
        }
    }
    names = st_tempname(cols(R))
    (void) st_addvar("double", names)
    st_store(., names, touse, R)
    st_local(outlocal, invtokens(names))
}
end
