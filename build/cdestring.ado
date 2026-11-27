*! version 1.4.0 20261127 github.com/mdroste/stata-ctools
*! original: Stata's destring (destring.ado 2.4.4, September 15, 2016), originally written by Nicholas J. Cox
*! original: see help cdestring for details

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

    * Both modes stage values and metadata in temporary variables. A plugin
    * error or Break then removes only those, without publishing partial
    * outputs. The variables are created after the plugin has parsed the
    * values and reported the storage type each one needs.
    unab original_order : _all
    local dest_varlist ""
    local i = 0
    foreach v of local varlist {
        local ++i
        tempvar tmp`i'
        local dest_varlist "`dest_varlist' `tmp`i''"
    }

    * Build option strings for C plugin
    local threads_code ""
    if `threads' > 0 {
        local threads_code "threads(`threads')"
    }

    * Keep user bytes out of Stata/plugin command syntax, including quotes.
    * The plugin gets every character it removes: ignore() and, with percent,
    * "%". It reports the ones a variable contained by their byte offsets in
    * this list.
    mata: st_local("__cds_ignore", st_local("ignore") + ///
            (st_local("percent") != "" & !strpos(st_local("ignore"), "%") ? "%" : "")); ///
        st_local("ignore_code", "ignorelen=" + strtrim(strofreal(strlen(st_local("__cds_ignore")), "%21.0f")))

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
        }
    }
    local plugin_nvars : word count `plugin_vars'

    * The plugin varlist holds the sources, in order: an index list in the
    * argument exceeds Stata's plugin-call argument length on wide data
    if `plugin_nvars' _ctools_strw `plugin_vars'

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

    * The plugin converts in two calls. parse reads the strings, keeps the
    * numbers and sets locals cdestring_failed (nonnumeric count),
    * cdestring_types (the storage type compress would choose),
    * cdestring_missing (missing values among the observations read) and
    * cdestring_seen (the ignore() characters found), one word per variable.
    * Each destination is then created with its final type and store fills
    * it, so Stata does not count, scan or recast the new variable. Break
    * waits until the plugin has released the numbers.
    local cdestring_failed ""
    local cdestring_types ""
    local cdestring_missing ""
    local cdestring_seen ""
    local cdestring_ignore_read ""
    local __time_create = 0
    nobreak {
        if `plugin_nvars' {
            capture noisily plugin call ctools_plugin `plugin_vars' `marker' `ifcond' `in', ///
                "cdestring parse `threads_code' nvars=`plugin_nvars' `force_code' `percent_code' `dpcomma_code' `float_code' `verbose_code' `select_code' `ignore_flags' `ignore_code'"
            local plugin_rc = _rc
            if `plugin_rc' exit `plugin_rc'
            if "`cdestring_ignore_read'" != "1" {
                capture plugin call ctools_plugin, "cdestring discard"
                di as error "cdestring: install the matching ctools plugin before converting data"
                exit 198
            }
            local n_converted = scalar(_cdestring_n_converted)
            local n_failed = scalar(_cdestring_n_failed)
            local nobs_read = scalar(_cdestring_nobs)
        }
        else {
            local n_converted = 0
            local n_failed = 0
            * (no plugin phase to time)
            foreach phase in parse load convert store total {
                scalar _cdestring_time_`phase' = 0
            }
        }

        * Create the destinations. Variables converted in Stata below start
        * as double (float under float) and are compressed afterwards.
        if `__do_timing' quietly _ctools_clock __ct_create0
        local plugin_publish
        local i = 0
        local j = 0
        capture noisily {
            foreach v of local varlist {
                local ++i
                if `: list v in stata_vars' {
                    local from_plugin_`i' 0
                    if "`float'" != "" quietly generate float `tmp`i'' = .
                    else quietly generate double `tmp`i'' = .
                }
                else {
                    local ++j
                    local from_plugin_`i' 1
                    local ptype_`i' : word `j' of `cdestring_types'
                    local pmiss_`i' : word `j' of `cdestring_missing'
                    local pseen_`i' : word `j' of `cdestring_seen'
                    local pseen_`i' : subinstr local pseen_`i' "," " ", all
                    if "`ptype_`i''" != "." {
                        quietly generate `ptype_`i'' `tmp`i'' = .
                        local plugin_publish `plugin_publish' `tmp`i''
                    }
                }
            }
        }
        local create_rc = _rc
        if `__do_timing' {
            quietly _ctools_clock __ct_now
            local __time_create = `__ct_now' - `__ct_create0'
        }
        if `create_rc' {
            if `plugin_nvars' capture plugin call ctools_plugin, "cdestring discard"
            exit `create_rc'
        }
        if "`plugin_publish'" != "" {
            capture noisily plugin call ctools_plugin `plugin_publish', "cdestring store `threads_code'"
            if _rc exit _rc
        }
        else if `plugin_nvars' {
            capture plugin call ctools_plugin, "cdestring discard"
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
            * (only a variable converted in Stata has a destination by now)
            if !`from_plugin_`i'' drop `dv'
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
                if `from_plugin_`i'' {
                    * The plugin listed the characters this variable held by
                    * where each first occurs in the list it was given.
                    mata: st_local("__cds_at", ///
                        strofreal(strpos(st_local("__cds_ignore"), st_local("t")) - 1))
                    local found : list __cds_at in pseen_`i'
                }
                else {
                    local find = cond(`asbytes', "strpos", "ustrpos")
                    if "`iftouse'" != "" quietly count if `find'(`v', `"`t'"') & `iftouse' `in'
                    else quietly count if `find'(`v', `"`t'"') `in'
                    local found = r(N) > 0
                }
                if `found' {
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
            * destring reports every missing value of the new variable and
            * compresses it to the smallest storage type. The plugin counted
            * the missing values among the observations it read (the others
            * are missing too) and the variable already has that type.
            if `from_plugin_`i'' {
                local nmiss_`i' = _N - `nobs_read' + `pmiss_`i''
                local type_`i' `ptype_`i''
            }
            else {
                quietly count if missing(`dv')
                local nmiss_`i' = r(N)
                quietly compress `dv'
                local type_`i' : type `dv'
            }
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
                * (order rewrites the whole dataset: skip it when the new
                * variables already stand where the old ones did)
                unab current_order : _all
                if "`current_order'" != "`original_order'" order `original_order'
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
                local target_order : list retokenize target_order
                unab current_order : _all
                if "`current_order'" != "`target_order'" order `target_order'
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
        capture local __plugin_time_total = scalar(_cdestring_time_total)
        if _rc != 0 local __plugin_time_total = 0
        local __plugin_call_overhead = `__time_plugin' - `__plugin_time_total' - `__time_create'

        di as text ""
        di as text "{hline 55}"
        di as text "cdestring timing breakdown:"
        di as text "{hline 55}"
        di as text "  C plugin internals:"
        di as text "    Argument parsing:       " as result %8.4f scalar(_cdestring_time_parse) " sec"
        di as text "    Data load:              " as result %8.4f scalar(_cdestring_time_load) " sec"
        di as text "    String conversion:      " as result %8.4f scalar(_cdestring_time_convert) " sec"
        di as text "    Store results:          " as result %8.4f scalar(_cdestring_time_store) " sec"
        di as text "  {hline 53}"
        di as text "    C plugin total:         " as result %8.4f scalar(_cdestring_time_total) " sec"
        di as text ""
        di as text "  Stata overhead:"
        di as text "    Pre-plugin setup:       " as result %8.4f `__time_preplugin' " sec"
        di as text "    Create variables:       " as result %8.4f `__time_create' " sec"
        di as text "    Plugin call overhead:   " as result %8.4f `__plugin_call_overhead' " sec"
        di as text "    Post-plugin cleanup:    " as result %8.4f `__time_postplugin' " sec"
        di as text "  {hline 53}"
        local __stata_overhead = `__time_preplugin' + `__time_create' + `__plugin_call_overhead' + `__time_postplugin'
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
        capture local __threads_max = scalar(_cdestring_threads_max)
        if _rc == 0 {
            capture local __openmp_enabled = scalar(_cdestring_openmp_enabled)
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
