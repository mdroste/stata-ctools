*! version 1.2.6 20261111 github.com/mdroste/stata-ctools

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

    * Mark sample
    marksample touse, novarlist

    * Count observations
    quietly count if `touse'
    local nobs_sample = r(N)
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

    * Create variable for keep flags (will be dropped at end)
    tempvar __keep
    quietly gen double `__keep' = 0

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
    if `"`if'"' == "" & `"`in'"' == "" local opts "`opts' skipif"
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

    * Call plugin with ALL variables (cencode pattern)
    * Pass keep_idx so C code knows which variable to write to
    plugin call ctools_plugin `allvars' `touse', ///
        "csample `threads_opt' `keep_idx' `nby' `by_indices' `opts'"

    * Drop observations where __keep == 0
    quietly drop if `touse' & `__keep' == 0

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
        di as text "    Data load:              " as result %8.4f _csample_time_load " sec"
        if `nby' > 0 {
            di as text "    Sort by-variables:      " as result %8.4f _csample_time_sort " sec"
        }
        di as text "    Sampling:               " as result %8.4f _csample_time_sample " sec"
        di as text "    Data store:             " as result %8.4f _csample_time_store " sec"
        di as text "  {hline 58}"
        di as text "    C plugin total:         " as result %8.4f _csample_time_total " sec"
        di as text "{hline 60}"
        di as text "    Wall clock total:       " as result %8.4f `__time_total' " sec"
        di as text "{hline 60}"
        di as text ""
        di as text "  Original obs: " as result `nobs_sample'
        di as text "  Kept:         " as result %8.0f _csample_kept
        di as text "  Dropped:      " as result %8.0f _csample_dropped
        di as text "  Groups:       " as result %8.0f _csample_ngroups
        di as text "  Threads:      " as result %8.0f _csample_threads
        di as text "{hline 60}"
    }

    * Summary message
    if "`verbose'" == "" {
        local kept = _csample_kept
        local dropped = _csample_dropped
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
