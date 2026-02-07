*! version 0.9.1 06Feb2026
*! cdecode: C-accelerated numeric to string decoding for Stata
*! Part of the ctools suite
*!
*! Replaces Stata's built-in decode command with a high-performance
*! C implementation.
*!
*! Syntax: cdecode varlist [if] [in], generate(newvarlist) | replace [maxlength(#) Verbose THReads(integer)]
*!
*! The cdecode command inherits the same syntax and functionality as Stata's
*! built-in decode command, producing identical results but with better
*! performance on large datasets.
*!
*! ctools extensions (not in Stata's decode):
*!   - varlist support: decode multiple numeric variables at once
*!   - replace option: replace original variables with decoded versions

* Mata function to parse label save file and write labels to output file
* This avoids Stata's command line length limits when there are many labels
capture mata: mata drop _cdecode_parse_labels_to_file()
mata:
void _cdecode_parse_labels_to_file(string scalar filename, string scalar lblname, string scalar outfile)
{
    real scalar fh_in, fh_out, val, lbl_len, max_label_len
    string scalar line, prefix, rest, lbl, lbl_escaped, entry
    string scalar open_cq, close_cq
    real scalar prefix_len, open_pos, close_pos, space_pos

    // Compound quote markers: `" and "'
    open_cq = char(96) + char(34)   // `"
    close_cq = char(34) + char(39)  // "'

    max_label_len = 0

    fh_in = fopen(filename, "r")
    if (fh_in < 0) {
        errprintf("cdecode: could not open label file\n")
        exit(error(601))
    }

    fh_out = fopen(outfile, "w")
    if (fh_out < 0) {
        fclose(fh_in)
        errprintf("cdecode: could not open output file\n")
        exit(error(601))
    }

    while ((line = fget(fh_in)) != J(0,0,"")) {
        // Line format: label define lblname VALUE `"TEXT"', modify
        prefix = "label define " + lblname + " "
        prefix_len = strlen(prefix)

        if (substr(line, 1, prefix_len) != prefix) {
            continue
        }

        rest = substr(line, prefix_len + 1, .)

        // Extract value (first token before space)
        space_pos = strpos(rest, " ")
        if (space_pos == 0) continue
        val = strtoreal(substr(rest, 1, space_pos - 1))
        rest = strtrim(substr(rest, space_pos + 1, .))

        // Find compound quotes `"..."'
        open_pos = strpos(rest, open_cq)
        close_pos = strpos(rest, close_cq)

        if (open_pos > 0 && close_pos > open_pos) {
            // Extract text between compound quotes
            lbl = substr(rest, open_pos + 2, close_pos - open_pos - 2)
        }
        else if (substr(rest, 1, 1) == char(34)) {
            // Simple quote format: "text", modify
            // Find closing quote before ", modify"
            close_pos = strpos(substr(rest, 2, .), char(34))
            if (close_pos > 0) {
                lbl = substr(rest, 2, close_pos - 1)
            }
            else {
                continue
            }
        }
        else {
            continue
        }

        // Escape special characters: \ -> \\, | -> \|, " -> \"
        lbl_escaped = subinstr(lbl, "\", "\\", .)
        lbl_escaped = subinstr(lbl_escaped, "|", "\|", .)
        lbl_escaped = subinstr(lbl_escaped, char(34), "\" + char(34), .)

        // Track max length
        lbl_len = strlen(lbl)
        if (lbl_len > max_label_len) {
            max_label_len = lbl_len
        }

        // Write entry to file (value|label on each line)
        entry = strofreal(val) + "|" + lbl_escaped
        fput(fh_out, entry)
    }

    fclose(fh_in)
    fclose(fh_out)

    // Return max label length to Stata
    st_local("max_label_len", strofreal(max_label_len))
}
end

program define cdecode, rclass
    version 14.1

    * Check observation limit (Stata plugin API limitation)
    if _N > 2147483647 {
        di as error "ctools does not support datasets exceeding 2^31 (2.147 billion) observations"
        di as error "This is a limitation of Stata's plugin API"
        exit 920
    }

    syntax varlist [if] [in], [Generate(string) replace MAXLength(integer 0) Verbose THReads(integer 0)]

    * =========================================================================
    * UPFRONT VALIDATION
    * =========================================================================

    * Check that either generate or replace is specified, but not both
    if "`generate'" == "" & "`replace'" == "" {
        di as error "cdecode: must specify either generate() or replace option"
        exit 198
    }
    if "`generate'" != "" & "`replace'" != "" {
        di as error "cdecode: cannot specify both generate() and replace options"
        exit 198
    }

    * Count input variables
    local n_vars : word count `varlist'

    * Handle replace vs generate
    local __do_replace = 0
    if "`replace'" != "" {
        local __do_replace = 1
    }
    else {
        * Check generate list has same number of variables
        local n_gen : word count `generate'
        if `n_gen' != `n_vars' {
            di as error "cdecode: generate() must specify `n_vars' variable(s) to match varlist"
            exit 198
        }
        * Check that none of the generate variables already exist
        forvalues i = 1/`n_gen' {
            local gvar : word `i' of `generate'
            capture confirm variable `gvar'
            if _rc == 0 {
                di as error "cdecode: variable `gvar' already exists"
                exit 110
            }
        }
    }

    * Check that all source variables are numeric and have value labels
    foreach v of local varlist {
        capture confirm numeric variable `v'
        if _rc != 0 {
            di as error "`v' is not a numeric variable"
            exit 108
        }
        local lblname : value label `v'
        if "`lblname'" == "" {
            di as error "cdecode: `v' has no value label attached"
            exit 182
        }
    }

    * =========================================================================
    * END VALIDATION
    * =========================================================================

    * Start timing
    local __do_timing = ("`verbose'" != "")
    if `__do_timing' {
        timer clear 90
        timer clear 91
        timer clear 92
        timer clear 93
        timer on 90
        timer on 91
    }

    * Mark sample
    marksample touse

    * Load the platform-appropriate ctools plugin
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
        if `__is_arm' == 0 {
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

    capture program ctools_plugin, plugin using("`__plugin'")
    if _rc != 0 & _rc != 110 & "`__plugin'" != "ctools.plugin" {
        capture program ctools_plugin, plugin using("ctools.plugin")
    }
    if _rc != 0 & _rc != 110 {
        di as error "cdecode: Could not load ctools plugin"
        exit 601
    }

    * Reset _rc
    capture confirm number 1

    * Build threads option
    local threads_code ""
    if `threads' > 0 {
        local threads_code "threads(`threads')"
    }

    * =========================================================================
    * Process each variable
    * =========================================================================

    forvalues i = 1/`n_vars' {
        local srcvar : word `i' of `varlist'

        * Determine destination variable name
        if `__do_replace' {
            tempvar __tempgen
            local destvar "`__tempgen'"
        }
        else {
            local destvar : word `i' of `generate'
        }

        * Get value label name attached to source variable
        local lblname : value label `srcvar'

        * Get variable index for source
        unab allvars : *
        local var_idx = 0
        local idx = 1
        foreach v of local allvars {
            if ("`v'" == "`srcvar'") {
                local var_idx = `idx'
                continue, break
            }
            local ++idx
        }

        if `var_idx' == 0 {
            di as error "cdecode: could not find variable `srcvar'"
            exit 111
        }

        * =====================================================================
        * Extract value labels and encode for C plugin
        * =====================================================================

        * Use label save to preserve trailing spaces in labels
        tempfile __lblfile
        quietly label save `lblname' using `__lblfile', replace

        * Create temp file for labels
        tempfile __labelsfile

        * Parse the label save file using Mata for reliable string handling
        mata: _cdecode_parse_labels_to_file("`__lblfile'", "`lblname'", "`__labelsfile'")

        * Determine string variable width
        if `maxlength' > 0 {
            local strwidth = `maxlength'
        }
        else {
            local strwidth = `max_label_len'
            if `strwidth' < 1 {
                local strwidth = 1
            }
        }

        * Cap at max string length
        if `strwidth' > 2045 {
            local strwidth = 2045
        }

        * Create the destination string variable
        quietly generate str`strwidth' `destvar' = ""

        * Build maxlen option
        local maxlen_code "maxlen=`strwidth'"

        * Call the C plugin
        if `__do_timing' {
            timer off 91
            timer on 92
        }
        plugin call ctools_plugin `srcvar' `destvar' `if' `in', ///
            `"cdecode `threads_code' `maxlen_code' labelsfile=`__labelsfile'"'
        if `__do_timing' {
            timer off 92
            timer on 93
        }

        local plugin_rc = _rc

        if `plugin_rc' != 0 {
            exit `plugin_rc'
        }

        * Handle replace option: drop original var and rename temp var
        if `__do_replace' {
            capture drop `srcvar'
            rename `destvar' `srcvar'
        }

        if `__do_timing' {
            timer off 93
            timer on 91
        }
    }

    if `__do_timing' {
        timer off 91
    }

    * Store rclass results
    return scalar N_vars = `n_vars'

    if `__do_timing' {
        timer off 90

        * Extract timer values
        quietly timer list 90
        local __time_total = r(t90)

        * Extract sub-timer values
        quietly timer list 91
        local __time_preplugin = r(t91)
        quietly timer list 92
        local __time_plugin = r(t92)
        quietly timer list 93
        local __time_postplugin = r(t93)

        * Calculate plugin call overhead
        capture local __plugin_time_total = _cdecode_time_total
        if _rc != 0 local __plugin_time_total = 0
        local __plugin_call_overhead = `__time_plugin' - `__plugin_time_total'

        di as text ""
        di as text "{hline 55}"
        di as text "cdecode timing breakdown:"
        di as text "{hline 55}"
        di as text "  C plugin internals:"
        di as text "    Argument parsing:       " as result %8.4f _cdecode_time_parse " sec"
        di as text "    Decode:                 " as result %8.4f _cdecode_time_decode " sec"
        di as text "  {hline 53}"
        di as text "    C plugin total:         " as result %8.4f _cdecode_time_total " sec"
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
        di as text "  Variables decoded:        " as result `n_vars'

        * Display thread diagnostics
        capture local __threads_max = _cdecode_threads_max
        if _rc == 0 {
            capture local __openmp_enabled = _cdecode_openmp_enabled
            if _rc != 0 local __openmp_enabled = 0
            di as text ""
            di as text "  Thread diagnostics:"
            di as text "    OpenMP enabled:         " as result %8.0f `__openmp_enabled'
            di as text "    Max threads available:  " as result %8.0f `__threads_max'
            di as text "{hline 55}"
        }

        * Clean up timing scalars
        capture scalar drop _cdecode_time_parse _cdecode_time_decode _cdecode_time_total
        capture scalar drop _cdecode_threads_max _cdecode_openmp_enabled
    }
end
