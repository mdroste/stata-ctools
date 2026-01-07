*! version 1.1.0 05Jan2026
*! cexport: C-accelerated delimited text export for Stata
*! Part of the ctools suite
*!
*! Description:
*!   High-performance replacement for export delimited using a C plugin
*!   with parallel data loading and chunked formatting.
*!
*! Syntax:
*!   cexport delimited [varlist] using filename [if] [in], [options]
*!
*! Options:
*!   delimiter(string)   - Field delimiter (default: ",")
*!   novarnames          - Don't write variable names as header row
*!   quote               - Quote all string fields
*!   noquoteif           - Don't automatically quote strings containing delimiters
*!   replace             - Overwrite existing file
*!   datafmt             - Use display formats for numeric variables (not yet implemented)
*!   verbose             - Display progress information

program define cexport, rclass
    version 14.0

    * Parse the subcommand
    gettoken subcmd 0 : 0, parse(" ,")

    if "`subcmd'" != "delimited" {
        di as error "cexport: unknown subcommand `subcmd'"
        di as error "Syntax: cexport delimited [varlist] using filename [, options]"
        exit 198
    }

    * Parse the rest with export delimited-style syntax
    syntax [varlist] using/ [if] [in], [Delimiter(string) NOVARNames QUOTE ///
        NOQUOTEif REPLACE DATAfmt Verbose TIMEit]

    * Mark sample
    marksample touse, novarlist

    * Validate file
    if "`replace'" == "" {
        capture confirm file `"`using'"'
        if !_rc {
            di as error "file `using' already exists"
            di as error "use the replace option to overwrite"
            exit 602
        }
    }

    * Set default delimiter
    if `"`delimiter'"' == "" {
        local delimiter ","
    }

    * Handle tab delimiter
    if `"`delimiter'"' == "tab" | `"`delimiter'"' == "\t" {
        local delimiter "	"
        local delim_display "tab"
    }
    else {
        local delim_display `"`delimiter'"'
    }

    * Validate delimiter
    if length(`"`delimiter'"') != 1 {
        di as error "cexport: delimiter must be a single character"
        exit 198
    }

    * If no varlist specified, use all variables
    if "`varlist'" == "" {
        qui ds
        local varlist `r(varlist)'
    }

    * Count variables and observations
    local nvars : word count `varlist'
    qui count `if' `in'
    local nobs = r(N)

    if `nvars' == 0 {
        di as error "cexport: no variables to export"
        exit 198
    }

    if `nobs' == 0 {
        di as error "cexport: no observations to export"
        exit 2000
    }

    * Load the platform-appropriate ctools plugin if not already loaded
    capture program list ctools_plugin
    if _rc != 0 {
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
            di as error "cexport: Could not load ctools plugin"
            exit 601
        }
    }

    * Display info if verbose
    if "`verbose'" != "" {
        di as text ""
        di as text "{hline 60}"
        di as text "cexport delimited: High-Performance CSV Export"
        di as text "{hline 60}"
        di as text "File:       " as result `"`using'"'
        di as text "Variables:  " as result `nvars'
        di as text "Obs:        " as result `nobs'
        di as text "Delimiter:  " as result "`delim_display'"
        di as text "Header:     " as result cond("`novarnames'" == "", "yes", "no")
        di as text "Quoting:    " as result cond("`quote'" != "", "all strings", cond("`noquoteif'" != "", "none", "as needed"))
        di as text "{hline 60}"
        di ""
    }

    * Build plugin arguments
    local opt_noheader = cond("`novarnames'" != "", "noheader", "")
    local opt_quote = cond("`quote'" != "", "quote", "")
    local opt_noquoteif = cond("`noquoteif'" != "", "noquoteif", "")
    local opt_verbose = cond("`verbose'" != "", "verbose", "")

    * Pass variable names to the plugin via global macro
    global _cexport_varnames `varlist'

    * Record start time
    timer clear 99
    timer on 99

    * Call the C plugin
    * Plugin expects: filename delimiter [options]
    capture noisily plugin call ctools_plugin `varlist' `if' `in', ///
        "cexport `using' `delimiter' `opt_noheader' `opt_quote' `opt_noquoteif' `opt_verbose'"

    local export_rc = _rc

    * Clean up global macro
    macro drop _cexport_varnames

    if `export_rc' {
        di as error "Error exporting data (rc=`export_rc')"
        exit `export_rc'
    }

    timer off 99
    quietly timer list 99
    local elapsed = r(t99)

    * Display summary
    di as text ""
    di as text "file " as result `"`using'"' as text " saved"
    di as text "(" as result %12.0fc `nobs' as text " observations written)"

    * Display timing if requested
    if "`timeit'" != "" | "`verbose'" != "" {
        di as text ""
        di as text "Timing breakdown:"
        di as text "  Load:  " as result %8.4f _cexport_time_load " sec"
        di as text "  Write: " as result %8.4f _cexport_time_write " sec"
        di as text "  Total: " as result %8.4f _cexport_time_total " sec"

        if `elapsed' > 0 {
            local rows_per_sec = `nobs' / `elapsed'
            di as text "  Speed: " as result %12.0fc `rows_per_sec' as text " rows/sec"
        }
    }

    timer clear

    * Return results
    return scalar N = `nobs'
    return scalar k = `nvars'
    return scalar time = `elapsed'
    return local filename `"`using'"'

end
