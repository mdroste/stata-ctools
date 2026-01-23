*! version 1.0.0 22Jan2026
*! cdecode: C-accelerated numeric to string decoding for Stata
*! Part of the ctools suite
*!
*! Replaces Stata's built-in decode command with a high-performance
*! C implementation.
*!
*! Syntax: cdecode varname [if] [in], generate(newvar) [maxlength(#) Verbose THReads(integer)]
*!
*! The cdecode command inherits the same syntax and functionality as Stata's
*! built-in decode command, producing identical results but with better
*! performance on large datasets.

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
    real scalar first_entry

    // Compound quote markers: `" and "'
    open_cq = char(96) + char(34)   // `"
    close_cq = char(34) + char(39)  // "'

    max_label_len = 0
    first_entry = 1

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
        if (first_entry) {
            fput(fh_out, entry)
            first_entry = 0
        }
        else {
            fput(fh_out, entry)
        }
    }

    fclose(fh_in)
    fclose(fh_out)

    // Return max label length to Stata
    st_local("max_label_len", strofreal(max_label_len))
}
end

program define cdecode
    version 14.0

    syntax varname(numeric) [if] [in], Generate(name) [MAXLength(integer 0) Verbose THReads(integer 0)]

    * =========================================================================
    * UPFRONT VALIDATION
    * =========================================================================

    * Check that source variable is numeric
    capture confirm numeric variable `varlist'
    if _rc != 0 {
        di as error "cdecode: `varlist' must be a numeric variable"
        exit 198
    }

    * Check that generate variable doesn't already exist
    capture confirm new variable `generate'
    if _rc != 0 {
        capture confirm variable `generate'
        if _rc == 0 {
            di as error "cdecode: variable `generate' already exists"
            exit 110
        }
        else {
            di as error "cdecode: invalid variable name `generate'"
            exit 198
        }
    }

    * Get value label name attached to source variable
    local lblname : value label `varlist'
    if "`lblname'" == "" {
        di as error "cdecode: `varlist' has no value label attached"
        exit 182
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

    * Get variable index
    unab allvars : *
    local var_idx = 0
    local idx = 1
    foreach v of local allvars {
        if ("`v'" == "`varlist'") {
            local var_idx = `idx'
            continue, break
        }
        local ++idx
    }

    if `var_idx' == 0 {
        di as error "cdecode: could not find variable `varlist'"
        exit 111
    }

    * =========================================================================
    * Extract value labels and encode for C plugin
    * =========================================================================

    * Use label save to preserve trailing spaces in labels
    * (The `: label` extended macro function strips trailing spaces)
    tempfile __lblfile
    quietly label save `lblname' using `__lblfile', replace

    * Build encoded label string: "value|label||value|label||..."
    * We need to escape | and \ in label text
    local labels_encoded ""
    local max_label_len = 0

    * Create temp file for labels (used when label string is too long)
    tempfile __labelsfile

    * Parse the label save file using Mata for reliable string handling
    * Format: label define lblname value `"text"', modify
    * Mata writes labels to temp file to avoid command line length limits
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
    quietly generate str`strwidth' `generate' = ""

    * Get index of new variable
    unab allvars : *
    local gen_idx = 0
    local idx = 1
    foreach v of local allvars {
        if ("`v'" == "`generate'") {
            local gen_idx = `idx'
            continue, break
        }
        local ++idx
    }

    * Build threads option
    local threads_code ""
    if `threads' > 0 {
        local threads_code "threads(`threads')"
    }

    * Build maxlen option
    local maxlen_code "maxlen=`strwidth'"

    if `__do_timing' {
        timer off 91
        timer on 92
    }

    * Call the C plugin (pass labels file path to avoid command line length limits)
    plugin call ctools_plugin `varlist' `generate' `if' `in', ///
        `"cdecode `threads_code' `maxlen_code' labelsfile=`__labelsfile'"'

    local plugin_rc = _rc

    if `__do_timing' {
        timer off 92
        timer on 93
    }

    * Get results
    local n_decoded = _cdecode_n_decoded
    local n_missing = _cdecode_n_missing
    local n_unlabeled = _cdecode_n_unlabeled

    if `__do_timing' {
        timer off 93
        timer off 90

        * Extract timer values
        quietly timer list 90
        local __time_total = r(t90)
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
        di as text "    Decode values:          " as result %8.4f _cdecode_time_decode " sec"
        di as text "  {hline 53}"
        di as text "    C plugin total:         " as result %8.4f _cdecode_time_total " sec"
        di as text "  {hline 53}"
        di as text "  Stata overhead:"
        di as text "    Pre-plugin setup:       " as result %8.4f `__time_preplugin' " sec"
        di as text "    Plugin call overhead:   " as result %8.4f `__plugin_call_overhead' " sec"
        di as text "    Post-plugin cleanup:    " as result %8.4f `__time_postplugin' " sec"
        di as text "{hline 55}"
        di as text "    Wall clock total:       " as result %8.4f `__time_total' " sec"
        di as text "{hline 55}"
        di as text ""
        di as text "  Values decoded:           " as result `n_decoded'
        di as text "  Missing values:           " as result `n_missing'
        if `n_unlabeled' > 0 {
            di as text "  Unlabeled values:         " as result `n_unlabeled'
        }
        di as text "  Max label length:         " as result `max_label_len'
        di as text "  String variable width:    " as result `strwidth'

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
    }

    if `plugin_rc' != 0 {
        exit `plugin_rc'
    }

    exit 0
end
