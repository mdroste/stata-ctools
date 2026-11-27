*! version 1.4.0 20261127 github.com/mdroste/stata-ctools
*! original: Stata's sample (sample.ado 2.4.0, January 20, 2021)
*! original: see help csample for details

program define csample, byable(onecall)
    version 14.1

    * Under version 6 or earlier, sample draws with its old algorithm, which
    * scripts rely on to reproduce old samples: sample itself does the draw
    if _caller() < 7 {
        local c = string(_caller())
        if _by() {
            di as error "sample may not be combined with by when run under version `c'"
            exit 190
        }
        _csample_version6 `c' `0'
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

    * sample # [, count] as in sample, where count makes # a number of
    * observations; count(#) is the older ctools form of the same thing.
    * pduplicates() sets how many uniforms sample draws to avoid ties; the
    * ctools draw (a partial shuffle) has no ties, so it is accepted unused.
    local cmdline `"`0'"'
    syntax [anything(name=percent)] [if] [in] [, *]
    local paren = ustrregexm(`"`options'"', "(^|\s)(c|co|cou|coun|count)\s*\(")
    local 0 `"`cmdline'"'
    if `paren' {
        syntax [anything(name=percent)] [if] [in], [Count(numlist integer max=1 >=0) ///
            BY(varlist) Verbose THReads(integer 0) PDUPlicates(real 1e-4)]
        if "`percent'" != "" {
            di as error "csample: cannot specify both # and count()"
            exit 198
        }
    }
    else {
        syntax [anything(name=percent)] [if] [in], [Count ///
            BY(varlist) Verbose THReads(integer 0) PDUPlicates(real 1e-4)]
        if "`percent'" == "" {
            di as error "csample: must specify # or count(#)"
            exit 198
        }
        * # is read as sample reads it (=/exp), with its errors
        _csample_exp `percent'
        if "`count'" != "" {
            local count `percent'
            local percent
        }
    }
    if !inrange(`threads', 0, 256) {
        di as error "threads() must be between 0 and 256"
        exit 198
    }
    if _N == 0 {
        foreach field in nobs kept dropped ngroups threads time_load time_sort time_sample time_store time_total {
            scalar _csample_`field' = 0
        }
    }
    * by varlist: is by(varlist), as in sample
    if _by() {
        if "`by'" != "" {
            di as error "by varlist: and by() option may not be combined"
            exit 190
        }
        local by `_byvars'
    }

    local __do_timing = ("`verbose'" != "")

    * Validate # as sample does, with its messages and return codes
    local use_percent = 0
    local use_count = 0
    local pct_value = 0

    if "`percent'" != "" {
        confirm number `percent'
        if `percent' < 0 | `percent' >= 100 {
            if `percent' == 100 exit
            di as error "{p 0 4}"
            di as error "`percent':"
            di as error ///
"# must be between 0 and 100 and is interpreted as a percent unless" _n ///
"option count is also specified, in which case # is interpreted as" _n ///
"the number of observations to be drawn" _n ///
"{p_end}"
            error 198
        }
        local use_percent = 1
        * The literal, as sample computes int(N*#/100+.5) from it
        local pct_value `percent'
    }
    else {
        confirm integer number `count'
        * sample prints this for a negative count and then keeps none
        if `count' < 0 {
            di as error "{p 0 4}"
            di as error "`count':"
            di as error ///
"# must be 0 or positive; with the count option, it represents the" _n ///
"# of observations to be drawn" _n ///
"{p_end}"
            local count 0
        }
        * A count of at least _N keeps everything, as in sample
        if `count' >= _N exit
        local use_count = 1
    }

    * Mark and count the sample. Without if and in it is the whole dataset:
    * the plugin is told so (skipif) and reads no marker.
    local whole = (`"`if'"' == "" & `"`in'"' == "")
    if `whole' {
        local touse
        local nobs_sample = _N
    }
    else {
        marksample touse, novarlist
        quietly count if `touse'
        local nobs_sample = r(N)
    }
    if `nobs_sample' == 0 {
        foreach field in nobs kept dropped ngroups time_load time_sort time_sample time_store time_total {
            scalar _csample_`field' = 0
        }
        scalar _csample_threads = 0
        exit 0
    }
    * sample leaves the data sorted by the by() variables
    local sort_by `by'

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

    * Count by-variables
    local nby = 0
    if "`by'" != "" {
        local nby : word count `by'
    }

    * strL by() values the plugin cannot read (binary, or over 2,045 bytes)
    * are grouped through stand-in ranks in their sort order
    if "`by'" != "" {
        local stubs
        foreach v of local by {
            tempvar stub
            local stubs `stubs' `stub'
        }
        _ctools_strlkey `by', generate(`stubs')
        local by `__ctools_keys'
    }

    * Create variable for keep flags (will be dropped at end). The plugin
    * writes the flag of every observation of the sample, so for the whole
    * dataset the variable needs no initial values.
    tempvar __keep
    if `whole' mata: (void) st_addvar("byte", st_local("__keep"), 1)
    else quietly gen double `__keep' = 0

    * Get ALL variables and find indices (cencode pattern). The sample
    * marker goes last, in place of an if: Stata answers an if one
    * observation at a time, while the plugin reads the marker in parallel.
    unab allvars : *
    local allvars : list allvars - touse
    local keep_idx = 0
    local idx = 1
    foreach v of local allvars {
        if ("`v'" == "`__keep'") {
            local keep_idx = `idx'
        }
        local ++idx
    }

    * Build by-variable indices
    local by_indices ""
    if "`by'" != "" {
        foreach bv of local by {
            local idx = 1
            foreach v of local allvars {
                if ("`v'" == "`bv'") {
                    local by_indices "`by_indices' `idx'"
                    continue, break
                }
                local ++idx
            }
        }
    }

    * Build options
    local opts ""
    if `use_percent' == 1 {
        local opts "percent=`pct_value'"
    }
    else if `count' == 0 {
        * Keeping no observation per group is a zero percent sample
        local opts "percent=0"
    }
    else {
        local count_s = strtrim(string(`count', "%20.0f"))
        local opts "count=`count_s'"
    }
    local opts "`opts' seedhi=`seed_hi' seedlo=`seed_lo'"
    if "`verbose'" != "" local opts "`opts' verbose"
    if `whole' local opts "`opts' skipif"
    else local opts "`opts' touse"
    * compact: the plugin moves a small sample to the top of the dataset
    * and the rest is dropped by range, where drop if evaluates every
    * observation and rewrites what is left. Not with e(sample) in memory,
    * whose hidden column the plugin cannot move with the observations.
    local esample : e(functions)
    if `whole' & !`: list posof "sample" in esample' {
        local opts "`opts' compact"
        _ctools_strw `allvars'
    }

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

    capture scalar drop _csample_compacted
    * Call plugin with ALL variables (cencode pattern)
    * Pass keep_idx so C code knows which variable to write to
    plugin call ctools_plugin `allvars' `touse', ///
        "csample `threads_opt' `keep_idx' `nby' `by_indices' `opts'"

    * Drop observations where __keep == 0
    local compacted = 0
    capture local compacted = scalar(_csample_compacted)
    if `whole' & `compacted' == 1 {
        * the kept observations are the first ones
        if scalar(_csample_kept) < _N quietly drop in `=scalar(_csample_kept) + 1'/l
    }
    else if `whole' quietly drop if `__keep' == 0
    else quietly drop if `touse' & `__keep' == 0

    * Clean up the keep variable
    drop `__keep'
    * As sample, leave the data sorted by the by() variables (stable: the
    * kept observations keep their order within each group)
    if "`sort_by'" != "" sort `sort_by', stable

    * Finish timing
    quietly _ctools_clock __ct_now
    local __ct_elapsed90 = `__ct_elapsed90' + (`__ct_now' - `__ct_start90')

    if `__do_timing' {
        * Elapsed seconds are held in caller-local timing state.
        local __time_total = `__ct_elapsed90'

        di as text ""
        di as text "{hline 60}"
        di as text "csample timing breakdown:"
        di as text "{hline 60}"
        di as text "  Stata pre-plugin:         " as result %8.4f `__time_preplugin' " sec"
        di as text "  {hline 58}"
        di as text "  C plugin phases:"
        di as text "    Data load:              " as result %8.4f scalar(_csample_time_load) " sec"
        if `nby' > 0 {
            di as text "    Sort by-variables:      " as result %8.4f scalar(_csample_time_sort) " sec"
        }
        di as text "    Sampling:               " as result %8.4f scalar(_csample_time_sample) " sec"
        di as text "    Data store:             " as result %8.4f scalar(_csample_time_store) " sec"
        di as text "  {hline 58}"
        di as text "    C plugin total:         " as result %8.4f scalar(_csample_time_total) " sec"
        di as text "{hline 60}"
        di as text "    Wall clock total:       " as result %8.4f `__time_total' " sec"
        di as text "{hline 60}"
        di as text ""
        di as text "  Original obs: " as result `nobs_sample'
        di as text "  Kept:         " as result %8.0f scalar(_csample_kept)
        di as text "  Dropped:      " as result %8.0f scalar(_csample_dropped)
        di as text "  Groups:       " as result %8.0f scalar(_csample_ngroups)
        di as text "  Threads:      " as result %8.0f scalar(_csample_threads)
        di as text "{hline 60}"
    }

    * Summary message
    if "`verbose'" == "" {
        local kept = scalar(_csample_kept)
        local dropped = scalar(_csample_dropped)
        di as text "(`dropped' " plural(`dropped', "observation") " deleted)"
    }
end

* sample under version # (6 or earlier), without the ctools-only options
program define _csample_version6
    gettoken c 0 : 0
    syntax [anything(everything)] [, Verbose THReads(integer 0) *]
    if `"`options'"' != "" local options `", `options'"'
    version `c', missing: sample `anything' `options'
end

* Check # as sample's syntax =/exp does (in its own scope: syntax resets the
* caller's if and in)
program define _csample_exp
    local 0 `"= `0'"'
    syntax =/exp
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

* _ctools_strw — set string width metadata for ctools flat buffer optimization
* 
* Usage:  _ctools_strw varlist
* 
* Sets local __ctools_strw in CALLER's scope containing comma-separated
* string variable widths (0 = numeric, >0 = str width e.g. 17 for str17).
* The ctools C plugin reads this local automatically via SF_macro_use
* to enable flat buffer string I/O (eliminates per-string strlen, memcpy,
* and arena atomic CAS overhead).
* 
* Example:
* _ctools_strw `allvars'
* plugin call ctools_plugin `allvars', "mycommand args"
* * C plugin automatically detects __ctools_strw local

program define _ctools_strw
    syntax varlist

    local __strw ""
    local __sep ""
    foreach v of local varlist {
        local __type : type `v'
        if substr("`__type'", 1, 3) == "str" {
            local __w = real(substr("`__type'", 4, .))
            if missing(`__w') local __w = 0
        }
        else {
            local __w = 0
        }
        local __strw "`__strw'`__sep'`__w'"
        local __sep ","
    }
    c_local __ctools_strw `__strw'
end
* ctools-helpers:end
