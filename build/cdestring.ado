*! version 1.2.6 20261111 github.com/mdroste/stata-ctools

program define cdestring
    version 14.1
    local cmdargs `"`0'"'

    * Check observation limit (Stata plugin API limitation)
    if _N > 2147483647 {
        di as error "ctools does not support datasets exceeding 2^31 (2.147 billion) observations"
        di as error "This is a limitation of Stata's plugin API"
        exit 920
    }

    * Parse syntax. As in destring, no varlist means all variables; numeric
    * ones are reported "already numeric" and skipped.
    syntax [varlist] [if] [in], [Generate(string) replace ///
        Ignore(string asis) force float percent dpcomma Verbose THReads(integer 0)]
    * The if is evaluated now, before any helper replaces the r() or e()
    * results it may name.
    local iftouse ""
    if `"`if'"' != "" {
        tempvar iftouse
        mark `iftouse' `if' `in'
        local if "if `iftouse'"
    }
    if !inrange(`threads', 0, 256) {
        di as error "threads() must be between 0 and 256"
        exit 198
    }

    * ignore("chars" [, aschars asbytes illegal]) is parsed as destring does:
    * a comma must be quoted, and _parse_ignore validates the suboptions.
    local ignore_flags ""
    if `"`ignore'"' != "" {
        if ustrpos(`"`ignore'"', ",") & !ustrpos(`"`ignore'"', `"""') {
            di as error `"ignore string must be enclosed with quotes to ignore commas or use options"'
            exit 198
        }
        _parse_ignore `ignore'
        local ignore `"`s(ignore)'"'
        local ignore_user `"`ignore'"'
        if !`s(aschars)' local ignore_flags "asbytes"
        * destring adds U+FFFD, which is also how Stata's Unicode functions
        * read each byte of an invalid sequence; the plugin removes those.
        if `s(illegal)' {
            local ignore `"`ignore'`=uchar(65533)'"'
            local ignore_flags "`ignore_flags' illegal"
        }
        sreturn clear
    }

    * =========================================================================
    * UPFRONT VALIDATION
    * =========================================================================

    * Checks and messages in destring's order and words
    if "`generate'" != "" & "`replace'" != "" {
        di as error "options generate and replace are mutually exclusive"
        exit 198
    }
    if "`generate'" == "" & "`replace'" == "" {
        di as error "must specify either generate or replace option"
        exit 198
    }

    * generate() is a newvarlist, so ranges such as n1-n3 work; its names are
    * checked before their number, as destring does
    if "`generate'" != "" {
        local ct1 : word count `varlist'
        local save_varlist `varlist'
        local save_if `"`if'"'
        local save_in `"`in'"'
        local 0 `"`generate'"'
        capture syntax newvarlist
        if _rc {
            di as error "generate() contains existing variable(s) and/or illegal variable name(s)"
            exit _rc
        }
        local generate `varlist'
        local varlist `save_varlist'
        local if `"`save_if'"'
        local in `"`save_in'"'
        if `ct1' != `: word count `generate'' {
            di as error "number of variables in varlist must equal"
            di as error "number of variables in generate(newvarlist)"
            exit 198
        }
    }
    if "`varlist'" == "" {
        di as error "varlist required"
        exit 100
    }

    * Messages go out per variable, in the order typed (msg_vars). destring
    * converts a repeated variable once and then finds it numeric; staging it
    * twice would make the group rename fail.
    local msg_vars `varlist'
    if "`replace'" != "" local varlist : list uniq varlist
    local genrep = cond("`replace'" != "", "replace", "generate")

    * Already-numeric variables are skipped (reported in the message loop)
    local string_vars ""
    local string_gen ""
    local numeric_vars ""
    local _var_idx = 0
    foreach v of local varlist {
        local ++_var_idx
        capture confirm string variable `v'
        if _rc != 0 {
            local numeric_vars `numeric_vars' `v'
        }
        else {
            local string_vars "`string_vars' `v'"
            if "`generate'" != "" {
                local gv : word `_var_idx' of `generate'
                local string_gen "`string_gen' `gv'"
            }
        }
    }

    * Update varlist and generate to only include string variables
    local varlist "`string_vars'"
    if "`generate'" != "" {
        local generate "`string_gen'"
    }
    local nvars : word count `varlist'

    * If no string variables remain, exit successfully (like destring)
    if `nvars' == 0 {
        foreach v of local msg_vars {
            di as text "`v' already numeric; no " as result "`genrep'"
        }
        exit 0
    }

    * =========================================================================
    * END VALIDATION
    * =========================================================================

    * Start timing
    local __do_timing = ("`verbose'" != "")
    if `__do_timing' {
        local __ct_elapsed90 = 0
        local __ct_elapsed91 = 0
        local __ct_elapsed92 = 0
        local __ct_elapsed93 = 0
        if "`__ct_elapsed90'" == "" local __ct_elapsed90 = 0
        quietly _ctools_clock __ct_start90
        if "`__ct_elapsed91'" == "" local __ct_elapsed91 = 0
        quietly _ctools_clock __ct_start91
    }

    * Load the platform-appropriate ctools plugin (cached after first call)
    _ctools_load
    * Stata scopes plugin registrations to the calling ado program.
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    capture confirm number 0

    * Both modes stage values and metadata. A plugin error or Break then
    * removes only temporary variables, without publishing partial outputs.
    unab original_order : _all
    local dest_varlist ""
    local i = 0
    foreach v of local varlist {
        local ++i
        tempvar tmp`i'
        if "`float'" != "" {
            quietly generate float `tmp`i'' = .
        }
        else {
            quietly generate double `tmp`i'' = .
        }
        local dest_varlist "`dest_varlist' `tmp`i''"
    }

    * Build option strings for C plugin
    local threads_code ""
    if `threads' > 0 {
        local threads_code "threads(`threads')"
    }

    local ignore_code ""
    if `"`ignore'"' != "" {
        * Escape special characters for passing to C
        * Escape backslashes first so literal \t differs from an actual tab.
        local ignore_escaped = subinstr(`"`ignore'"', "\", "\\", .)
        local ignore_escaped = subinstr(`"`ignore_escaped'"', " ", "\s", .)
        local ignore_escaped = subinstr(`"`ignore_escaped'"', char(9), "\t", .)
        local ignore_escaped = subinstr(`"`ignore_escaped'"', char(10), "\n", .)
        local ignore_escaped = subinstr(`"`ignore_escaped'"', char(13), "\r", .)
        local ignore_code "ignore=`ignore_escaped'"
    }

    local force_code ""
    if "`force'" != "" {
        local force_code "force"
    }

    local percent_code ""
    if "`percent'" != "" {
        local percent_code "percent"
    }

    local dpcomma_code ""
    if "`dpcomma'" != "" {
        local dpcomma_code "dpcomma"
    }

    * The plugin needs float for the float range and percent rounding
    local float_code ""
    if "`float'" != "" {
        local float_code "float"
    }

    local verbose_code ""
    if "`verbose'" != "" {
        local verbose_code "verbose"
    }

    if `__do_timing' {
        quietly _ctools_clock __ct_now
        local __ct_elapsed91 = `__ct_elapsed91' + (`__ct_now' - `__ct_start91')
        if "`__ct_elapsed92'" == "" local __ct_elapsed92 = 0
        quietly _ctools_clock __ct_start92
    }

    * strL variables holding values the plugin cannot read (binary, or text
    * over 2,045 bytes) are converted in Stata below; the plugin takes the rest
    _ctools_strlkey `varlist', check
    local stata_vars `__ctools_strl_bad'
    * Paired by position: generate() may convert one source twice
    local plugin_vars
    local plugin_dest
    local stata_dest
    local stata_src
    local i = 0
    foreach v of local varlist {
        local ++i
        if `: list v in stata_vars' {
            local stata_src `stata_src' `v'
            local stata_dest `stata_dest' `: word `i' of `dest_varlist''
        }
        else {
            local plugin_vars `plugin_vars' `v'
            local plugin_dest `plugin_dest' `: word `i' of `dest_varlist''
        }
    }
    local plugin_nvars : word count `plugin_vars'

    * The plugin varlist holds the sources, then the destinations (seqpairs):
    * an index list in the argument exceeds Stata's plugin-call argument
    * length on wide data
    if `plugin_nvars' _ctools_strw `plugin_vars' `plugin_dest'

    * Stata answers an if one observation at a time. Without one the plugin
    * skips the check; the if marked after syntax is passed last, as a
    * sample marker the plugin reads in parallel.
    local select_code "skipif"
    local marker ""
    local ifcond ""
    if `"`if'"' != "" {
        local marker "`iftouse'"
        local select_code "touse"
    }

    * ignore= goes last so its characters cannot be read as option keywords.
    * The plugin sets local cdestring_failed (nonnumeric count per variable).
    local cdestring_failed ""
    if `plugin_nvars' {
        capture noisily plugin call ctools_plugin `plugin_vars' `plugin_dest' `marker' `ifcond' `in', ///
            "cdestring `threads_code' seqpairs nvars=`plugin_nvars' `force_code' `percent_code' `dpcomma_code' `float_code' `verbose_code' `select_code' `ignore_flags' `ignore_code'"
        local plugin_rc = _rc
        if `plugin_rc' exit `plugin_rc'
        local n_converted = _cdestring_n_converted
        local n_failed = _cdestring_n_failed
    }
    else {
        local n_converted = 0
        local n_failed = 0
        * (no plugin phase to time)
        foreach phase in parse load convert store total {
            scalar _cdestring_time_`phase' = 0
        }
    }

    * The other strL variables, through destring's own steps: remove the
    * ignore() characters, then "%" (percent), then invalid UTF-8 (illegal),
    * trim, note the missing values, apply dpcomma, cut leading zeros to one
    * and take real(). Values real() leaves missing that were not missing
    * text are nonnumeric; with percent, any "%" divides the variable by 100.
    local stata_failed
    if "`stata_vars'" != "" {
        local chars `"`ignore_user'"'
        if "`percent'" != "" & !ustrpos(`"`chars'"', "%") local chars `"`chars'%"'
        if strpos("`ignore_flags'", "illegal") local chars `"`chars'`=uchar(65533)'"'
        local asbytes = strpos("`ignore_flags'", "asbytes") > 0
        local nchars = cond(`asbytes', length(`"`chars'"'), ustrlen(`"`chars'"'))
        local sample = cond("`iftouse'" != "", "if `iftouse'", "")
        tempvar text ismissing
        foreach v of local stata_src {
            gettoken dv stata_dest : stata_dest
            * Observations outside the sample keep "", a missing value
            quietly generate strL `text' = `v' `sample' `in'
            quietly count if strpos(`text', "%")
            local haspct = r(N) > 0 & "`percent'" != ""
            forvalues m = 1/`nchars' {
                if `asbytes' {
                    local t = substr(`"`chars'"', `m', 1)
                    quietly replace `text' = subinstr(`text', `"`t'"', "", .)
                }
                else {
                    local t = usubstr(`"`chars'"', `m', 1)
                    quietly replace `text' = usubinstr(`text', `"`t'"', "", .)
                }
            }
            quietly replace `text' = ustrtrim(`text')
            quietly generate byte `ismissing' = `text' == "" | `text' == "." | ///
                (strlen(`text') == 2 & inrange(`text', ".a", ".z"))
            if "`dpcomma'" != "" {
                quietly replace `text' = subinstr(`text', ".", "X", .) if !`ismissing'
                quietly replace `text' = subinstr(`text', ",", ".", 1) ///
                    if !ustrregexm(`text', "(,[a-z])|(^,$)")
            }
            quietly replace `text' = ustrregexrf(`text', "^0+", "0")
            quietly replace `dv' = real(`text')
            if `haspct' quietly replace `dv' = `dv' / 100
            quietly count if missing(`dv') & !`ismissing'
            local nbad = r(N)
            quietly count `sample' `in'
            local n_converted = `n_converted' + r(N) - `nbad'
            local n_failed = `n_failed' + `nbad'
            local stata_failed `stata_failed' `nbad'
            drop `text' `ismissing'
        }
    }

    * Nonnumeric counts in the order of varlist
    local failed_plugin `cdestring_failed'
    local cdestring_failed
    foreach v of local varlist {
        if `: list v in stata_vars' {
            gettoken nbad stata_failed : stata_failed
            local cdestring_failed `cdestring_failed' `nbad'
        }
        else {
            gettoken nbad failed_plugin : failed_plugin
            local cdestring_failed `cdestring_failed' `nbad'
        }
    }

    if `__do_timing' {
        quietly _ctools_clock __ct_now
        local __ct_elapsed92 = `__ct_elapsed92' + (`__ct_now' - `__ct_start92')
        if "`__ct_elapsed93'" == "" local __ct_elapsed93 = 0
        quietly _ctools_clock __ct_start93
    }

    * Like destring, a variable with any nonnumeric value (after removing
    * ignore() characters) is left untouched unless force is specified: no
    * generate, no replace (the plugin did not store it). Others proceed.
    * destring records the ignore() characters it removed in char destring
    local ignore_chars `"`ignore'"'
    if "`percent'" != "" & !ustrpos(`"`ignore_chars'"', "%") local ignore_chars `"`ignore_chars'%"'
    * As destring does, asbytes lists bytes and illegal is shown by name.
    local asbytes = strpos("`ignore_flags'", "asbytes") > 0
    local aswhat = cond(`asbytes', "byte", "character")
    local n_ignore = cond(`asbytes', length(`"`ignore_chars'"'), ustrlen(`"`ignore_chars'"'))
    local publish
    local final_names
    local backups
    local i = 1
    foreach v of local varlist {
        local dv : word `i' of `dest_varlist'
        local nbad : word `i' of `cdestring_failed'
        if "`nbad'" == "" local nbad 0
        local nbad_`i' `nbad'
        if `nbad' > 0 & "`force'" == "" {
            drop `dv'
            local done_`i' 0
        }
        else {
            local done_`i' 1
            local vlabel : variable label `v'
            label variable `dv' `"`vlabel'"'
            local removed ""
            local nremoved = 0
            local seen_chars ""
            forvalues m = 1/`n_ignore' {
                if `asbytes' local t = substr(`"`ignore_chars'"', `m', 1)
                else local t = usubstr(`"`ignore_chars'"', `m', 1)
                if strpos(`"`seen_chars'"', `"`t'"') continue
                local seen_chars `"`seen_chars'`t'"'
                local find = cond(`asbytes', "strpos", "ustrpos")
                if "`iftouse'" != "" quietly count if `find'(`v', `"`t'"') & `iftouse' `in'
                else quietly count if `find'(`v', `"`t'"') `in'
                if r(N) {
                    if `"`t'"' == uchar(65533) & strpos("`ignore_flags'", "illegal") {
                        local removed `"`removed' illegal Unicode"'
                    }
                    else if `"`t'"' == " " local removed `"`removed' space"'
                    else local removed `"`removed' `t'"'
                    local ++nremoved
                }
            }
            * destring shows the bytes outside ASCII as \xhh (printablebytes())
            if `asbytes' local removed = printablebytes(`"`removed'"')
            local removed_`i' `"`removed'"'
            local nremoved_`i' `nremoved'
            local removed_text = plural(`nremoved', cond(`asbytes', "Byte", "Character")) + " removed " + ///
                plural(`nremoved', "was", "were") + `":`removed'"'
            char `dv'[destring] `"`removed_text'"'
            char `dv'[destring_cmd] `"cdestring `cmdargs'"'
            char `v'[destring]
            char `v'[destring_cmd]
            * destring copies the characteristics for generate() (through
            * macros) and moves them with char rename for replace (below)
            if "`replace'" == "" {
                local characteristics : char `v'[]
                foreach ch of local characteristics {
                    local value : char `v'[`ch']
                    char `dv'[`ch'] `"`value'"'
                }
            }
            * destring reports every missing value of the new variable
            quietly count if missing(`dv')
            local nmiss_`i' = r(N)
            * Compress to the optimal storage type (mimics destring)
            quietly compress `dv'
            local type_`i' : type `dv'
            local publish `publish' `dv'
            if "`replace'" != "" {
                tempvar backup`i'
                local backups `backups' `backup`i''
                local final_names `final_names' `v'
            }
            else {
                local final : word `i' of `generate'
                local final_names `final_names' `final'
                local converted_sources `converted_sources' `v'
            }
        }
        local ++i
    }

    * destring's messages, one variable at a time in the order typed: a
    * repeated variable that was replaced is then already numeric
    local seen
    foreach v of local msg_vars {
        if `: list v in numeric_vars' | ("`replace'" != "" & `: list v in seen' & `: list v in final_names') {
            di as text "`v' already numeric; no " as result "`genrep'"
            continue
        }
        local seen `seen' `v'
        if "`replace'" != "" {
            local i : list posof "`v'" in varlist
            local finalprint ""
        }
        else {
            local i : word count `seen'
            local finalprint "`: word `i' of `generate'' "
        }
        if !`done_`i'' {
            if `"`ignore_chars'"' != "" {
                di as text "`v': contains `aswhat's not specified in " ///
                    as result "ignore()" as text "; no " as result "`genrep'"
            }
            else {
                di as text "`v': contains nonnumeric `aswhat's; no " as result "`genrep'"
            }
            continue
        }
        if `nbad_`i'' {
            di as text "`v': contains nonnumeric `aswhat's" _continue
            if `"`ignore_chars'"' != "" di " not specified in " as result "ignore()" _continue
            if `"`removed_`i''"' == "" {
                di "; `finalprint'" as result "`genrep'd " as text "as " as result "`type_`i''" _continue
            }
            di
        }
        if `"`removed_`i''"' != "" {
            di as text "`v': " plural(`nremoved_`i'', "`aswhat'") as result `"`removed_`i''"' ///
                as text " removed; `finalprint'" as result "`genrep'd " as text "as " as result "`type_`i''"
        }
        else if !`nbad_`i'' {
            di as text "`v': all characters numeric; `finalprint'" as result "`genrep'd " ///
                as text "as " as result "`type_`i''"
        }
        if `nmiss_`i'' {
            di as text "(`nmiss_`i'' missing " plural(`nmiss_`i'', "value") " generated)"
        }
    }

    * Keep originals until the complete name swap succeeds. Group rename
    * validates its whole mapping; nobreak prevents interruption mid-commit.
    if "`publish'" != "" {
        nobreak {
            if "`replace'" != "" rename (`final_names') (`backups')
            capture rename (`publish') (`final_names')
            local publish_rc = _rc
            if `publish_rc' {
                if "`replace'" != "" rename (`backups') (`final_names')
                exit `publish_rc'
            }
            if "`replace'" != "" {
                foreach backup of local backups {
                    gettoken final final_names : final_names
                    char rename `backup' `final'
                }
                drop `backups'
                order `original_order'
            }
            else {
                * destring places each generate() output right after its
                * source, so a later output of the same source comes first
                local target_order ""
                local repeated : list dups converted_sources
                foreach v of local original_order {
                    local target_order `target_order' `v'
                    local pos : list posof "`v'" in converted_sources
                    * (blocks: the word macro is invalid for pos 0, and an
                    * inline if expands it before testing)
                    if `pos' & !`: list v in repeated' {
                        local target_order `target_order' `: word `pos' of `final_names''
                    }
                    else if `pos' {
                        local outputs
                        local j = 0
                        foreach src of local converted_sources {
                            local ++j
                            if "`src'" == "`v'" local outputs `: word `j' of `final_names'' `outputs'
                        }
                        local target_order `target_order' `outputs'
                    }
                }
                order `target_order'
            }
        }
    }

    if `__do_timing' {
        quietly _ctools_clock __ct_now
        local __ct_elapsed93 = `__ct_elapsed93' + (`__ct_now' - `__ct_start93')
        quietly _ctools_clock __ct_now
        local __ct_elapsed90 = `__ct_elapsed90' + (`__ct_now' - `__ct_start90')

        * Extract timer values
        * Elapsed seconds are held in caller-local timing state.
        local __time_total = `__ct_elapsed90'
        * Elapsed seconds are held in caller-local timing state.
        local __time_preplugin = `__ct_elapsed91'
        * Elapsed seconds are held in caller-local timing state.
        local __time_plugin = `__ct_elapsed92'
        * Elapsed seconds are held in caller-local timing state.
        local __time_postplugin = `__ct_elapsed93'

        * Calculate plugin call overhead
        capture local __plugin_time_total = _cdestring_time_total
        if _rc != 0 local __plugin_time_total = 0
        local __plugin_call_overhead = `__time_plugin' - `__plugin_time_total'

        di as text ""
        di as text "{hline 55}"
        di as text "cdestring timing breakdown:"
        di as text "{hline 55}"
        di as text "  C plugin internals:"
        di as text "    Argument parsing:       " as result %8.4f _cdestring_time_parse " sec"
        di as text "    Data load:              " as result %8.4f _cdestring_time_load " sec"
        di as text "    String conversion:      " as result %8.4f _cdestring_time_convert " sec"
        di as text "    Store results:          " as result %8.4f _cdestring_time_store " sec"
        di as text "  {hline 53}"
        di as text "    C plugin total:         " as result %8.4f _cdestring_time_total " sec"
        di as text ""
        di as text "  Stata overhead:"
        di as text "    Pre-plugin setup:       " as result %8.4f `__time_preplugin' " sec"
        di as text "    Plugin call overhead:   " as result %8.4f `__plugin_call_overhead' " sec"
        di as text "    Post-plugin cleanup:    " as result %8.4f `__time_postplugin' " sec"
        di as text "  {hline 53}"
        local __stata_overhead = `__time_preplugin' + `__plugin_call_overhead' + `__time_postplugin'
        di as text "    Stata overhead total:   " as result %8.4f `__stata_overhead' " sec"
        di as text "{hline 55}"
        di as text "    Wall clock total:       " as result %8.4f `__time_total' " sec"
        di as text "{hline 55}"
        di as text ""
        di as text "  Numeric or missing:       " as result `n_converted'
        if `n_failed' > 0 {
            di as text "  Nonnumeric:               " as result `n_failed'
        }

        * Display thread diagnostics
        capture local __threads_max = _cdestring_threads_max
        if _rc == 0 {
            capture local __openmp_enabled = _cdestring_openmp_enabled
            if _rc != 0 local __openmp_enabled = 0
            di as text ""
            di as text "  Thread diagnostics:"
            di as text "    OpenMP enabled:         " as result %8.0f `__openmp_enabled'
            di as text "    Max threads available:  " as result %8.0f `__threads_max'
            di as text "{hline 55}"
        }
    }

    * Ensure clean exit
    exit 0
end
