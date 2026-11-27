*! version 1.3.2 20261124 github.com/mdroste/stata-ctools
*! original: Stata's bsample (bsample.ado 2.0.8, January 17, 2025)
*! original: see help cbsample for details

program define cbsample
    version 14.1
    * Legacy draws use bsample_7; parse the ctools-only options under 14.1.
    if _caller() < 8 {
        _cbsample_version7 `0'
        exit
    }

    * Check observation limit (Stata plugin API limitation)
    if _N > 2147483647 {
        di as error "ctools does not support datasets exceeding 2^31 (2.147 billion) observations"
        di as error "This is a limitation of Stata's plugin API"
        exit 920
    }

    local __ct_elapsed90 = 0
    if "`__ct_elapsed90'" == "" local __ct_elapsed90 = 0
    quietly _ctools_clock __ct_start90

    * Parse syntax to match bsample exactly:
    * bsample [exp] [if] [in], [strata(varlist) cluster(varlist) idcluster(newvar) weight(varname)]
    syntax [anything(name=n)] [if] [in], [STRata(varlist) CLuster(varlist) IDcluster(string) Weight(varname) Verbose THReads(integer 0)]
    local __wgtvar `weight'
    local weight
    _ctools_cleanup cbsample
    if !inrange(`threads', 0, 256) {
        di as error "threads() must be between 0 and 256"
        exit 198
    }

    local __do_timing = ("`verbose'" != "")

    * With no rows, native bsample distinguishes its fast no-option path
    * from explicit count/if/strata/cluster syntax. Let it validate those
    * empty-data cases; no random draws or C transfers are required.
    if _N == 0 {
        local empty_opts
        if "`strata'" != "" local empty_opts `empty_opts' strata(`strata')
        if "`cluster'" != "" local empty_opts `empty_opts' cluster(`cluster')
        if "`__wgtvar'" != "" local empty_opts `empty_opts' weight(`__wgtvar')
        if `"`idcluster'"' != "" local empty_opts `empty_opts' idcluster(`idcluster')
        quietly bsample `n' `if' `in', `empty_opts'
        foreach result in nobs n_selected total_weight nclusters nstrata {
            scalar _cbsample_`result' = 0
        }
        scalar _cbsample_threads = 1
        exit
    }

    local n_exp `"`n'"'

    * Validate weight variable exists if specified (matches bsample behavior;
    * with idcluster() bsample reports the idcluster() error first)
    if "`__wgtvar'" != "" & `"`idcluster'"' == "" {
        confirm numeric variable `__wgtvar'
    }

    * A stratified draw with no eligible strata has zero weight. An empty
    * unstratified selection in nonempty data is an error instead.
    marksample touse, novarlist
    markout `touse' `strata', strok
    * bsample leaves the caller's r() in place.
    tempname __rhold
    _return hold `__rhold'
    quietly count if `touse'
    local nobs_sample = r(N)
    _return restore `__rhold'
    * idcluster(), checked and an existing variable of that name dropped
    * before the draw, as bsample does. bsample has sorted the data by
    * sample, strata and cluster when it reports these errors.
    if `"`idcluster'"' != "" {
        local id_error
        if "`cluster'" == "" local id_error "option idcluster() requires option cluster()"
        else if "`__wgtvar'" != "" local id_error "options idcluster() and weight() may not be combined"
        else {
            capture confirm new variable `idcluster'
            if _rc {
                capture confirm variable `idcluster'
                if _rc local id_error confirm
            }
        }
        if "`id_error'" != "" {
            if `nobs_sample' < _N {
                tempvar presort
                quietly generate byte `presort' = cond(`touse', 1, .)
                sort `presort' `strata' `cluster', stable
                drop `presort'
            }
            else if "`strata'`cluster'" != "" sort `strata' `cluster', stable
            if "`id_error'" == "confirm" confirm variable `idcluster'
            di as error "`id_error'"
            exit 198
        }
        capture confirm new variable `idcluster'
        if _rc drop `idcluster'
    }
    * The resample size is an expression, evaluated as bsample does: on the
    * sample, within each stratum under strata() (so _N counts the stratum)
    * Sizes that differ across strata are stored by stratum in a variable.
    tempvar nsampvar
    if `"`n_exp'"' == "" | `nobs_sample' == 0 local n = 0
    else _cbsample_size n : `"`n_exp'"' `touse' "`strata'" `nsampvar' "`cluster'"
    if `nobs_sample' == 0 {
        if "`strata'" == "" {
            di as error "cbsample: no observations in sampling sample"
            exit 198
        }
        if "`__wgtvar'" != "" quietly replace `__wgtvar' = 0
        else quietly drop if 1
        if `"`idcluster'"' != "" {
            quietly generate byte `idcluster' = .
            label variable `idcluster' "Bootstrap sample cluster id"
            sort `strata'
        }
        foreach result in nobs n_selected total_weight nclusters nstrata {
            scalar _cbsample_`result' = 0
        }
        scalar _cbsample_threads = 1
        exit
    }

    * Load plugin
    _ctools_load
    * Stata scopes plugin registrations to the calling ado program.
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    capture confirm number 0

    * Generate seed from Stata's RNG (respects -set seed-). Pass two exact
    * 31-bit halves: a single ~62-bit value loses digits in a macro.
    local seed_hi = string(floor(runiform() * 2147483647), "%10.0f")
    local seed_lo = string(floor(runiform() * 2147483647), "%10.0f")

    * Count cluster and strata variables
    local ncluster = 0
    if "`cluster'" != "" {
        local ncluster : word count `cluster'
    }

    local nstrata = 0
    if "`strata'" != "" {
        local nstrata : word count `strata'
    }

    * strL strata() and cluster() values the plugin cannot read (binary, or
    * over 2,045 bytes) are grouped through stand-in ranks in their sort
    * order; Stata's own steps below keep the variables themselves
    local plugin_cluster `cluster'
    local plugin_strata `strata'
    foreach role in cluster strata {
        if "``role''" != "" {
            local stubs
            foreach v of local `role' {
                tempvar stub
                local stubs `stubs' `stub'
            }
            _ctools_strlkey ``role'', generate(`stubs')
            local plugin_`role' `__ctools_keys'
        }
    }

    * Create temporary frequency weight variable
    tempvar __cbsample_weight__
    quietly gen double `__cbsample_weight__' = 0
    * idcluster(): the draw also writes each sampled observation's sorted
    * position and its cluster and stratum ordinals
    if `"`idcluster'"' != "" {
        tempvar __cbsample_pos__ __cbsample_cl__ __cbsample_st__ __cbsample_key__
        foreach v in pos cl st {
            quietly gen double `__cbsample_`v'__' = .
        }
    }

    * Get ALL variables and find indices. The sample marker goes last, in
    * place of an if: Stata answers an if one observation at a time, while
    * the plugin reads the marker in parallel.
    unab allvars : *
    local allvars : list allvars - touse
    local freq_idx = 0
    local size_idx = 0
    local id_opts
    local idx = 1
    foreach v of local allvars {
        if ("`v'" == "`__cbsample_weight__'") {
            local freq_idx = `idx'
        }
        if "`n'" == "var" & "`v'" == "`nsampvar'" local size_idx = `idx'
        if `"`idcluster'"' != "" {
            if "`v'" == "`__cbsample_pos__'" local id_opts `id_opts' posidx=`idx'
            if "`v'" == "`__cbsample_cl__'" local id_opts `id_opts' clidx=`idx'
            if "`v'" == "`__cbsample_st__'" local id_opts `id_opts' stidx=`idx'
        }
        local ++idx
    }
    * bsample numbers the cluster copies across strata under version 14 or
    * later, within each stratum otherwise
    if `"`idcluster'"' != "" & _caller() > 13 local id_opts `id_opts' iduniq

    * Build cluster variable indices
    local cluster_indices ""
    if "`cluster'" != "" {
        foreach cv of local plugin_cluster {
            local idx = 1
            foreach v of local allvars {
                if ("`v'" == "`cv'") {
                    local cluster_indices "`cluster_indices' `idx'"
                    continue, break
                }
                local ++idx
            }
        }
    }

    * Build strata variable indices
    local strata_indices ""
    if "`strata'" != "" {
        foreach sv of local plugin_strata {
            local idx = 1
            foreach v of local allvars {
                if ("`v'" == "`sv'") {
                    local strata_indices "`strata_indices' `idx'"
                    continue, break
                }
                local ++idx
            }
        }
    }

    * Build options
    if "`n'" == "var" local n_s 0
    else local n_s = strtrim(string(`n', "%20.0f"))
    local opts "n=`n_s' seedhi=`seed_hi' seedlo=`seed_lo' `id_opts'"
    if `size_idx' local opts "`opts' sizeidx=`size_idx'"
    if "`verbose'" != "" local opts "`opts' verbose"
    if `nobs_sample' == _N local opts "`opts' skipif"
    else local opts "`opts' touse"

    local threads_opt ""
    if `threads' > 0 {
        local threads_opt "threads(`threads')"
    }

    * Record pre-plugin time
    quietly _ctools_clock __ct_now
    local __ct_elapsed90 = `__ct_elapsed90' + (`__ct_now' - `__ct_start90')
    * Elapsed seconds are held in caller-local timing state.
    local __time_preplugin = `__ct_elapsed90'
    if "`__ct_elapsed90'" == "" local __ct_elapsed90 = 0
    quietly _ctools_clock __ct_start90

    * Call plugin with ALL variables
    plugin call ctools_plugin `allvars' `touse', ///
        "cbsample `threads_opt' `freq_idx' `ncluster' `nstrata' `cluster_indices' `strata_indices' `opts'"

    * Record post-plugin time
    quietly _ctools_clock __ct_now
    local __ct_elapsed90 = `__ct_elapsed90' + (`__ct_now' - `__ct_start90')
    * Elapsed seconds are held in caller-local timing state.
    local __time_postplugin_start = `__ct_elapsed90'
    if "`__ct_elapsed90'" == "" local __ct_elapsed90 = 0
    quietly _ctools_clock __ct_start90

    * Handle output based on weight option (match bsample behavior exactly)
    if "`__wgtvar'" != "" {
        * weight() specified: copy weights to user's variable, keep all obs
        quietly replace `__wgtvar' = `__cbsample_weight__'
    }
    else {
        * No weight(): expand data to match bsample behavior
        if `"`idcluster'"' != "" {
            if scalar(_cbsample_total_weight) > 2147483647 {
                di as error "cbsample: expanded data exceed the plugin observation limit"
                exit 920
            }
            local idmax = scalar(_cbsample_idmax)
            local idtype = cond(`idmax' <= 100, "byte", cond(`idmax' <= 32740, "int", "long"))
            * (the temporary key first: if either variable cannot be added,
            * no idcluster() variable is left behind)
            quietly generate double `__cbsample_key__' = .
            quietly generate `idtype' `idcluster' = .
            format `idcluster' %8.0g
        }
        * No copy of the data is kept, as in bsample: an error from here on
        * leaves the observations dropped or expanded so far. Break waits
        * until the sample is complete.
        nobreak {
            * Drop observations with zero weight first
            quietly drop if `__cbsample_weight__' == 0
            * Expand based on weight (creates duplicates for weight > 1)
            quietly expand `__cbsample_weight__'
            if `"`idcluster'"' != "" {
                local id_pass "seedhi=`seed_hi' seedlo=`seed_lo'"
                if _caller() > 13 local id_pass "`id_pass' iduniq"
                if scalar(_cbsample_idrandom) local id_pass "`id_pass' idrandom"
                plugin call ctools_plugin `__cbsample_pos__' `__cbsample_cl__' `__cbsample_st__' ///
                    `idcluster' `__cbsample_key__', "cbsample idpass `id_pass'"
                label variable `idcluster' "Bootstrap sample cluster id"
                * bsample's row order: strata, id, then the observation order
                sort `strata' `__cbsample_key__'
            }
        }
    }

    * Finish timing
    quietly _ctools_clock __ct_now
    local __ct_elapsed90 = `__ct_elapsed90' + (`__ct_now' - `__ct_start90')

    if `__do_timing' {
        * Elapsed seconds are held in caller-local timing state.
        local __time_total = `__ct_elapsed90'
        local __time_postplugin = `__time_total' - `__time_postplugin_start'

        di as text ""
        di as text "{hline 60}"
        di as text "cbsample timing breakdown:"
        di as text "{hline 60}"
        di as text "  Stata pre-plugin:         " as result %8.4f `__time_preplugin' " sec"
        di as text "  {hline 58}"
        di as text "  C plugin phases:"
        di as text "    Data load:              " as result %8.4f scalar(_cbsample_time_load) " sec"
        if `ncluster' > 0 | `nstrata' > 0 {
            di as text "    Sort groups:            " as result %8.4f scalar(_cbsample_time_sort) " sec"
        }
        di as text "    Sampling:               " as result %8.4f scalar(_cbsample_time_sample) " sec"
        di as text "    Data store:             " as result %8.4f scalar(_cbsample_time_store) " sec"
        di as text "  {hline 58}"
        di as text "    C plugin total:         " as result %8.4f scalar(_cbsample_time_total) " sec"
        di as text "  {hline 58}"
        di as text "  Stata post-plugin:        " as result %8.4f `__time_postplugin' " sec"
        if "`__wgtvar'" == "" {
            di as text "    (drop zeros + expand)"
        }
        di as text "{hline 60}"
        di as text "    Wall clock total:       " as result %8.4f `__time_total' " sec"
        di as text "{hline 60}"
        di as text ""
        di as text "  Original obs: " as result `nobs_sample'
        di as text "  Selected:     " as result %8.0f scalar(_cbsample_n_selected)
        di as text "  Total weight: " as result %8.0f scalar(_cbsample_total_weight)
        di as text "  Clusters:     " as result %8.0f scalar(_cbsample_nclusters)
        di as text "  Strata:       " as result %8.0f scalar(_cbsample_nstrata)
        di as text "  Threads:      " as result %8.0f scalar(_cbsample_threads)
        di as text "{hline 60}"
    }
end

* bsample_7, without the ctools-only options
program define _cbsample_version7
    syntax [anything(everything)] [, Verbose THReads(integer 0) *]
    if `"`options'"' != "" local options `", `options'"'
    bsample_7 `anything' `options'
end

* The resample size bsample takes from [exp]: constant (within each stratum
* under strata()), an integer, and greater than zero. Stata evaluates _N in
* the sample (strata: in the stratum). Sizes that differ across strata are
* stored by observation in the new variable sizevar and macname is "var".
program define _cbsample_size
    args macname colon exp touse strata sizevar cluster
    capture confirm number `exp'
    if !_rc {
        local values `exp'
        local constant 1
    }
    else if "`strata'" == "" {
        * bsample evaluates subscripts after moving eligible rows first.
        * Preserve restores the caller's observation order on success or error.
        preserve
        tempvar nsamp eligible
        quietly generate byte `eligible' = !`touse'
        sort `eligible' `cluster', stable
        quietly generate double `nsamp' = (`exp') if `touse'
        mata: __cb_v = st_data(., st_local("nsamp"), st_local("touse")); ///
            st_local("constant", strofreal(all(__cb_v :== __cb_v[1]))); ///
            st_local("values", strofreal(__cb_v[1], "%21.0g"))
        mata: mata drop __cb_v
        restore
    }
    else {
        * bsample evaluates the expression by stratum, the data sorted by
        * (strata, cluster) with ties in their order
        quietly generate double `sizevar' = .
        preserve
        tempvar nsamp first obs
        quietly generate double `obs' = _n
        quietly keep if `touse'
        sort `strata' `cluster', stable
        quietly by `strata': generate double `nsamp' = (`exp')
        quietly by `strata': generate byte `first' = _n == 1
        quietly by `strata': replace `first' = . if `nsamp' != `nsamp'[1]
        quietly count if missing(`first')
        local constant = r(N) == 0
        mata: __cb_v = uniqrows(st_data(., st_local("nsamp"), st_local("first"))); ///
            st_local("values", invtokens(strofreal(__cb_v', "%21.0g")))
        if `constant' & `: word count `values'' > 1 {
            mata: __cb_s = st_data(., (st_local("obs"), st_local("nsamp")))
        }
        mata: mata drop __cb_v
        restore
        if !`constant' {
            capture drop `sizevar'
            di as error "expression: `exp': is not constant within strata"
            exit 198
        }
    }
    if !`constant' {
        di as error "expression: `exp': is not constant"
        exit 198
    }
    * bsample checks every size is an integer, then that every size is
    * greater than zero
    foreach value of local values {
        if !(abs(`value' - int(`value')) < 1e-7) {
            capture mata: mata drop __cb_s
            capture drop `sizevar'
            di as error "expression: `exp': is not an integer"
            exit 198
        }
    }
    foreach value of local values {
        if !(`value' > 0) {
            capture mata: mata drop __cb_s
            capture drop `sizevar'
            di as error "resample size must be greater than zero"
            exit 498
        }
    }
    if `: word count `values'' > 1 {
        mata: st_store(__cb_s[., 1], st_local("sizevar"), round(__cb_s[., 2]))
        mata: mata drop __cb_s
        c_local `macname' var
        exit
    }
    capture drop `sizevar'
    local value = round(`values')
    c_local `macname' `value'
end

* ctools-helpers:begin (generated by scripts/sync_release.py from src/ado; edit those files)
* ctools: drop the scalars and matrices a command and its plugin leave behind
* _ctools_cleanup <command> [<command> ...] drops every scalar and matrix
* named __<command>_* or _<command>_*.
program define _ctools_cleanup
    version 14.1
    foreach command of local 0 {
        foreach pattern in __`command'_* _`command'_* {
            local found : all scalars "`pattern'"
            if "`found'" != "" capture scalar drop `found'
            local found : all matrices "`pattern'"
            if "`found'" != "" capture matrix drop `found'
        }
    }
end

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

* _ctools_strlkey — numeric stand-ins for strL keys the plugin cannot read
* 
* Usage:  _ctools_strlkey varlist, generate(namelist) [missing]
* _ctools_strlkey varlist, check
* 
* The plugin reads strL text of at most 2,045 bytes and no binary values.
* Each strL variable of varlist that holds another value anywhere in the
* data gets a stand-in: the corresponding name of generate() becomes a long
* holding the rank of the value among the variable's distinct values, in
* Stata's sort order. Grouping by the ranks, ordering by them and counting
* their distinct values give what the strL values give. With missing, ""
* (a missing string) gets a missing rank, so that commands that leave out
* missing keys still leave it out; otherwise "" ranks first, as it sorts.
* 
* Sets local __ctools_keys in the caller: varlist with each such variable
* replaced by its stand-in. Variables the plugin can read stay as they are.
* The order of the observations, the sort marker and r() are left alone.
* check only sets local __ctools_strl_bad in the caller: the variables of
* varlist that would get a stand-in.

program define _ctools_strlkey
    version 14.1
    syntax varlist, [GENerate(namelist) MISSing CHECK]
    if "`check'" == "" & `: word count `generate'' != `: word count `varlist'' exit 198

    local strl
    foreach v of local varlist {
        if "`: type `v''" == "strL" local strl `strl' `v'
    }
    local strl : list uniq strl
    if "`strl'" == "" | _N == 0 {
        c_local __ctools_keys `varlist'
        c_local __ctools_strl_bad
        exit
    }

    * strlcheck reads lengths and binary flags only. If the plugin cannot
    * answer, every strL key gets a stand-in, which is exact anyway.
    local bad `strl'
    capture _ctools_load
    if !_rc {
        capture program ctools_plugin, plugin using("`__ctools_plugin'")
        if !_rc | _rc == 110 {
            local ctools_strl_bad ""
            capture plugin call ctools_plugin `strl', "strlcheck"
            if !_rc & `: word count `ctools_strl_bad'' == `: word count `strl'' {
                local bad
                local i = 0
                foreach v of local strl {
                    local ++i
                    if `: word `i' of `ctools_strl_bad'' local bad `bad' `v'
                }
            }
        }
    }

    if "`check'" != "" {
        c_local __ctools_strl_bad `bad'
        exit
    }

    local keys
    local made
    local i = 0
    foreach v of local varlist {
        local ++i
        if !`: list v in bad' {
            local keys `keys' `v'
            continue
        }
        * A variable named twice shares one stand-in
        local at : list posof "`v'" in made
        if `at' {
            local keys `keys' `: word `at' of `madeid''
            continue
        }
        local id : word `i' of `generate'
        _ctools_strlkey_rank `v' `id' `missing'
        local made `made' `v'
        local madeid `madeid' `id'
        local keys `keys' `id'
    }
    c_local __ctools_keys `keys'
end

* Rank the values of v in Stata's sort order (sortpreserve puts the
* observations, and the sort marker, back; on Break as well).
program define _ctools_strlkey_rank, sortpreserve
    args v id missing
    sort `v'
    quietly by `v': generate long `id' = _n == 1
    quietly replace `id' = sum(`id')
    if "`missing'" != "" quietly replace `id' = . if `v' == ""
end
* ctools-helpers:end
