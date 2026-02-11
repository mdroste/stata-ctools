*! version 1.0.2 9feb2026 github.com/mdroste/stata-ctools

program define cexport, rclass
    version 14.1

    * Check observation limit (Stata plugin API limitation)
    if _N > 2147483647 {
        di as error "ctools does not support datasets exceeding 2^31 (2.147 billion) observations"
        di as error "This is a limitation of Stata's plugin API"
        exit 920
    }

    * Parse the subcommand
    gettoken subcmd 0 : 0, parse(" ,")

    if "`subcmd'" == "excel" {
        * Dispatch to Excel export handler
        cexport_excel `0'
        return add
        exit
    }
    else if "`subcmd'" != "delimited" {
        di as error "cexport: unknown subcommand `subcmd'"
        di as error "Syntax: cexport delimited [varlist] [using] filename [, options]"
        di as error "        cexport excel [varlist] [using] filename [, options]"
        exit 198
    }

    * Check that data is loaded (match export delimited rc=102)
    if c(k) == 0 {
        di as error "no variables defined"
        exit 102
    }

    * Parse the rest with export delimited-style syntax
    * Support both "using filename" and just "filename" (like export delimited does)
    syntax [anything] [using/] [if] [in], [Delimiter(string) NOVARNames QUOTE ///
        NOQUOTEif REPLACE DATAfmt DATEString(string) NOLabel Verbose ///
        MMAP NOFSYNC DIRECT PREFAULT CRLF NOPARALLEL THReads(integer 0)]

    * Handle filename - can come from using/ or as last positional argument
    if `"`using'"' == "" & `"`anything'"' != "" {
        * No 'using' keyword: extract filename (token with dot) from anything
        local __tokens `"`anything'"'
        local __filename ""
        local __varlist ""
        while `"`__tokens'"' != "" {
            gettoken __tok __tokens : __tokens
            if strpos(`"`__tok'"', ".") > 0 {
                local __filename `"`__tok'"'
            }
            else {
                local __varlist `"`__varlist' `__tok'"'
            }
        }
        if `"`__filename'"' == "" {
            di as error "cexport delimited: filename required"
            di as error "Syntax: cexport delimited [varlist] [using] filename [, options]"
            exit 198
        }
        local using `"`__filename'"'
        local anything `"`__varlist'"'
    }
    if `"`using'"' == "" {
        di as error "cexport delimited: filename required"
        di as error "Syntax: cexport delimited [varlist] [using] filename [, options]"
        exit 198
    }
    * Resolve varlist from anything.
    * Must happen BEFORE marksample, which creates a touse tempvar
    * that ds would include in the variable list.
    if strtrim(`"`anything'"') != "" {
        unab varlist : `anything'
    }
    else {
        qui ds
        local varlist `r(varlist)'
    }

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

    * Handle tab delimiter - pass "tab" keyword to plugin (not actual char)
    local plugin_delim `"`delimiter'"'
    if `"`delimiter'"' == "tab" | `"`delimiter'"' == "\t" {
        local plugin_delim "tab"
        local delim_display "tab"
    }
    else {
        local delim_display `"`delimiter'"'
        * Validate delimiter (tab is handled specially)
        if length(`"`delimiter'"') != 1 {
            di as error "cexport: delimiter must be a single character"
            exit 198
        }
    }

    * Handle value labels: decode labeled variables unless nolabel specified
    local export_varlist ""
    local decoded_vars ""
    local original_names ""
    if "`nolabel'" == "" {
        * Default: export value labels
        foreach var of local varlist {
            * Check if variable has a value label attached
            local vallbl : value label `var'
            if "`vallbl'" != "" {
                * Has value label - decode to temporary string variable
                tempvar decoded_`var'
                qui decode `var', gen(`decoded_`var'')
                local export_varlist `export_varlist' `decoded_`var''
                local decoded_vars `decoded_vars' `decoded_`var''
                local original_names `original_names' `var'
            }
            else {
                local export_varlist `export_varlist' `var'
            }
        }
    }
    else {
        * nolabel specified: export raw numeric values
        local export_varlist `varlist'
    }

    * Handle datafmt and datestring(): format date/time variables as strings
    * This is done after value label handling (which may have created temp vars)
    if "`datafmt'" != "" | `"`datestring'"' != "" {
        local new_export_varlist ""
        local formatted_vars ""
        local export_idx = 0
        foreach var of local export_varlist {
            local export_idx = `export_idx' + 1
            * Get original variable name for format lookup
            * (export_varlist may contain temp vars from value label decoding)
            local origvar : word `export_idx' of `varlist'

            * Check if this is a date/time variable by looking at its format
            local fmt : format `origvar'
            local is_date_fmt = 0
            if substr("`fmt'", 1, 2) == "%t" | substr("`fmt'", 1, 3) == "%-t" {
                local is_date_fmt = 1
            }

            if `is_date_fmt' {
                * This is a date/time variable - convert to formatted string
                tempvar fmtvar_`export_idx'

                * Use datestring format if specified, otherwise use variable's display format
                * Note: missing values should become empty strings (like Stata's behavior)
                if `"`datestring'"' != "" {
                    qui gen str `fmtvar_`export_idx'' = cond(missing(`origvar'), "", string(`origvar', `"`datestring'"'))
                }
                else {
                    qui gen str `fmtvar_`export_idx'' = cond(missing(`origvar'), "", string(`origvar', "`fmt'"))
                }

                local new_export_varlist `new_export_varlist' `fmtvar_`export_idx''
                local formatted_vars `formatted_vars' `fmtvar_`export_idx''
            }
            else {
                * Not a date variable - keep as is
                local new_export_varlist `new_export_varlist' `var'
            }
        }
        local export_varlist `new_export_varlist'
    }

    * Count variables and observations
    local nvars : word count `varlist'
    qui count `if' `in'
    local nobs = r(N)

    if `nvars' == 0 {
        di as error "cexport: no variables to export"
        exit 198
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

    * Build plugin arguments
    local opt_noheader = cond("`novarnames'" != "", "noheader", "")
    local opt_quote = cond("`quote'" != "", "quote", "")
    local opt_noquoteif = cond("`noquoteif'" != "", "noquoteif", "")
    local opt_verbose = cond("`verbose'" != "", "verbose", "")

    * Performance options
    local opt_mmap = cond("`mmap'" != "", "mmap", "")
    local opt_nofsync = cond("`nofsync'" != "", "nofsync", "")
    local opt_direct = cond("`direct'" != "", "direct", "")
    local opt_prefault = cond("`prefault'" != "", "prefault", "")
    local opt_crlf = cond("`crlf'" != "", "crlf", "")
    local opt_noparallel = cond("`noparallel'" != "", "noparallel", "")

    * Build threads option string
    local threads_code ""
    if `threads' > 0 {
        local threads_code "threads(`threads')"
    }

    * Pass variable names to the plugin via global macro
    global CEXPORT_VARNAMES `varlist'

    * Pass variable storage types to the plugin
    * Types: 1=byte, 2=int, 3=long, 4=float, 5=double, 0=string
    local vartypes ""
    foreach var of local export_varlist {
        local vtype : type `var'
        if "`vtype'" == "byte" {
            local vartypes `vartypes' 1
        }
        else if "`vtype'" == "int" {
            local vartypes `vartypes' 2
        }
        else if "`vtype'" == "long" {
            local vartypes `vartypes' 3
        }
        else if "`vtype'" == "float" {
            local vartypes `vartypes' 4
        }
        else if "`vtype'" == "double" {
            local vartypes `vartypes' 5
        }
        else {
            * String type (str#, strL)
            local vartypes `vartypes' 0
        }
    }
    global CEXPORT_VARTYPES `vartypes'

    * Record start time
    timer clear 99
    timer on 99

    * Set string width metadata for flat buffer optimization
    _ctools_strw `export_varlist'

    * Call the C plugin
    * Plugin expects: filename delimiter [options]
    * Use export_varlist (may contain decoded temp vars for value labels)
    capture noisily plugin call ctools_plugin `export_varlist' `if' `in', ///
        "cexport `threads_code' `using' `plugin_delim' `opt_noheader' `opt_quote' `opt_noquoteif' `opt_verbose' `opt_mmap' `opt_nofsync' `opt_direct' `opt_prefault' `opt_crlf' `opt_noparallel'"

    local export_rc = _rc

    * Clean up global macros
    macro drop CEXPORT_VARNAMES
    macro drop CEXPORT_VARTYPES

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
    if "`verbose'" != "" {
        capture local __time_load = _cexport_time_load
        if _rc != 0 local __time_load = .
        capture local __time_format = _cexport_time_format
        if _rc != 0 local __time_format = .
        capture local __time_write = _cexport_time_write
        if _rc != 0 local __time_write = .
        capture local __time_total = _cexport_time_total
        if _rc != 0 local __time_total = .

        if `__time_total' != . {
            local plugin_overhead = `elapsed' - `__time_total'
            if `plugin_overhead' < 0 local plugin_overhead = 0
        }
        else {
            local plugin_overhead = .
        }
        di as text ""
        di as text "{hline 55}"
        di as text "cexport timing breakdown:"
        di as text "{hline 55}"
        di as text "  C plugin internals:"
        di as text "    Data load:              " as result %8.4f `__time_load' " sec"
        di as text "    Format:                 " as result %8.4f `__time_format' " sec"
        di as text "    Write to file:          " as result %8.4f `__time_write' " sec"
        di as text "  {hline 53}"
        di as text "    C plugin total:         " as result %8.4f `__time_total' " sec"
        di as text ""
        di as text "  Stata overhead:"
        di as text "    Plugin call overhead:   " as result %8.4f `plugin_overhead' " sec"
        di as text "  {hline 53}"
        di as text "    Stata overhead total:   " as result %8.4f `plugin_overhead' " sec"
        di as text "{hline 55}"
        di as text "    Wall clock total:       " as result %8.4f `elapsed' " sec"
        di as text "{hline 55}"
        if `elapsed' > 0 {
            local rows_per_sec = `nobs' / `elapsed'
            di as text "    Throughput:             " as result %12.0fc `rows_per_sec' as text " rows/sec"
        }

        * Display thread diagnostics
        capture local __threads_max = _cexport_threads_max
        if _rc == 0 {
            capture local __openmp_enabled = _cexport_openmp_enabled
            if _rc != 0 local __openmp_enabled = 0
            di as text ""
            di as text "  Thread diagnostics:"
            di as text "    OpenMP enabled:         " as result %8.0f `__openmp_enabled'
            di as text "    Max threads available:  " as result %8.0f `__threads_max'
            di as text "{hline 55}"
        }
    }

    * Return results
    return scalar N = `nobs'
    return scalar k = `nvars'
    return scalar time = `elapsed'
    return local filename `"`using'"'

    timer clear 99
end

/*******************************************************************************
 * cexport_excel: Export to Excel (.xlsx) files
 *
 * Syntax: cexport excel [varlist] [using] filename [if] [in], [options]
 *
 * Options:
 *   sheet(name)    - Worksheet name (default: "Sheet1")
 *   firstrow(variables|nonames) - First row handling
 *   replace        - Overwrite existing file
 *   nolabel        - Export numeric values instead of value labels
 *   verbose        - Display timing information
 ******************************************************************************/
program define cexport_excel, rclass
    version 14.1

    * Check observation limit (Stata plugin API limitation)
    if _N > 2147483647 {
        di as error "ctools does not support datasets exceeding 2^31 (2.147 billion) observations"
        di as error "This is a limitation of Stata's plugin API"
        exit 920
    }

    * Support both "using filename" and just "filename" (like export excel does)
    syntax [anything] [using/] [if] [in], [SHEET(string) FIRSTrow(string) ///
        REPLACE DATAfmt DATEString(string) NOLabel Verbose ///
        CELL(string) MISSING(string) KEEPCELLfmt]

    * Handle filename - can come from using/ or as last positional argument
    if `"`using'"' == "" & `"`anything'"' != "" {
        local __tokens `"`anything'"'
        local __filename ""
        local __varlist ""
        while `"`__tokens'"' != "" {
            gettoken __tok __tokens : __tokens
            if strpos(`"`__tok'"', ".") > 0 {
                local __filename `"`__tok'"'
            }
            else {
                local __varlist `"`__varlist' `__tok'"'
            }
        }
        if `"`__filename'"' == "" {
            di as error "cexport excel: filename required"
            di as error "Syntax: cexport excel [varlist] [using] filename [, options]"
            exit 198
        }
        local using `"`__filename'"'
        local anything `"`__varlist'"'
    }
    if `"`using'"' == "" {
        di as error "cexport excel: filename required"
        di as error "Syntax: cexport excel [varlist] [using] filename [, options]"
        exit 198
    }
    * Resolve varlist from anything.
    * Must happen BEFORE marksample, which creates a touse tempvar
    * that ds would include in the variable list.
    if strtrim(`"`anything'"') != "" {
        unab varlist : `anything'
    }
    else {
        qui ds
        local varlist `r(varlist)'
    }

    * Mark sample
    marksample touse, novarlist

    * Validate cell option format (e.g., A1, B5, AA10)
    if `"`cell'"' != "" {
        local cell_upper = upper(`"`cell'"')
        * Check format: 1-3 letters followed by 1 or more digits
        if !regexm("`cell_upper'", "^[A-Z]+[0-9]+$") {
            di as error "cexport excel: cell() must be a valid cell reference (e.g., A1, B5)"
            exit 198
        }
        local cell `"`cell_upper'"'
    }

    * Validate file extension
    local ext = substr(lower(`"`using'"'), -5, 5)
    if "`ext'" != ".xlsx" {
        di as error "cexport excel: file must have .xlsx extension"
        exit 198
    }

    * Validate file
    if "`replace'" == "" {
        capture confirm file `"`using'"'
        if !_rc {
            di as error "file `using' already exists"
            di as error "use the replace option to overwrite"
            exit 602
        }
    }

    * Parse firstrow option
    local opt_firstrow = "firstrow"
    if "`firstrow'" != "" {
        if "`firstrow'" == "nonames" {
            local opt_firstrow = "nofirstrow"
        }
        else if "`firstrow'" != "variables" {
            di as error "cexport excel: firstrow() must be 'variables' or 'nonames'"
            exit 198
        }
    }

    * Handle value labels
    local export_varlist ""
    local decoded_vars ""
    if "`nolabel'" == "" {
        foreach var of local varlist {
            local vallbl : value label `var'
            if "`vallbl'" != "" {
                tempvar decoded_`var'
                qui decode `var', gen(`decoded_`var'')
                local export_varlist `export_varlist' `decoded_`var''
                local decoded_vars `decoded_vars' `decoded_`var''
            }
            else {
                local export_varlist `export_varlist' `var'
            }
        }
    }
    else {
        local export_varlist `varlist'
    }

    * Handle datafmt and datestring(): format date/time variables as strings
    if "`datafmt'" != "" | `"`datestring'"' != "" {
        local new_export_varlist ""
        local formatted_vars ""
        local export_idx = 0
        foreach var of local export_varlist {
            local export_idx = `export_idx' + 1
            * Get original variable name for format lookup
            local origvar : word `export_idx' of `varlist'

            * Check if this is a date/time variable by looking at its format
            local fmt : format `origvar'
            local is_date_fmt = 0
            if substr("`fmt'", 1, 2) == "%t" | substr("`fmt'", 1, 3) == "%-t" {
                local is_date_fmt = 1
            }

            if `is_date_fmt' {
                * This is a date/time variable - convert to formatted string
                tempvar fmtvar_`export_idx'

                * Use datestring format if specified, otherwise use variable's display format
                * Note: missing values should become empty strings (like Stata's behavior)
                if `"`datestring'"' != "" {
                    qui gen str `fmtvar_`export_idx'' = cond(missing(`origvar'), "", string(`origvar', `"`datestring'"'))
                }
                else {
                    qui gen str `fmtvar_`export_idx'' = cond(missing(`origvar'), "", string(`origvar', "`fmt'"))
                }

                local new_export_varlist `new_export_varlist' `fmtvar_`export_idx''
                local formatted_vars `formatted_vars' `fmtvar_`export_idx''
            }
            else {
                * Not a date variable - keep as is
                local new_export_varlist `new_export_varlist' `var'
            }
        }
        local export_varlist `new_export_varlist'
    }
    else {
        * No datafmt/datestring: convert date values to Excel serial numbers
        * so they are stored as native Excel dates (not strings)
        local new_export_varlist ""
        local date_converted_vars ""
        local export_idx = 0
        local has_date_cols = 0
        foreach var of local export_varlist {
            local export_idx = `export_idx' + 1
            local origvar : word `export_idx' of `varlist'
            local fmt : format `origvar'
            local is_date_fmt = 0
            if substr("`fmt'", 1, 2) == "%t" | substr("`fmt'", 1, 3) == "%-t" {
                local is_date_fmt = 1
            }

            if `is_date_fmt' {
                local has_date_cols = 1
                * Determine date type: strip % and optional -
                local fmtsub = substr("`fmt'", 2, .)
                if substr("`fmtsub'", 1, 1) == "-" local fmtsub = substr("`fmtsub'", 2, .)
                local datetype = substr("`fmtsub'", 2, 1)

                tempvar xldate_`export_idx'
                if "`datetype'" == "d" {
                    qui gen double `xldate_`export_idx'' = `origvar' + 21916 if !missing(`origvar')
                }
                else if "`datetype'" == "m" {
                    qui gen double `xldate_`export_idx'' = dofm(`origvar') + 21916 if !missing(`origvar')
                }
                else if "`datetype'" == "q" {
                    qui gen double `xldate_`export_idx'' = dofq(`origvar') + 21916 if !missing(`origvar')
                }
                else if "`datetype'" == "w" {
                    qui gen double `xldate_`export_idx'' = dofw(`origvar') + 21916 if !missing(`origvar')
                }
                else if "`datetype'" == "h" {
                    qui gen double `xldate_`export_idx'' = dofh(`origvar') + 21916 if !missing(`origvar')
                }
                else if "`datetype'" == "y" {
                    qui gen double `xldate_`export_idx'' = dofy(`origvar') + 21916 if !missing(`origvar')
                }
                else if "`datetype'" == "c" | "`datetype'" == "C" {
                    if "`datetype'" == "C" {
                        qui gen double `xldate_`export_idx'' = dofC(`origvar') + 21916 if !missing(`origvar')
                    }
                    else {
                        qui gen double `xldate_`export_idx'' = dofc(`origvar') + 21916 if !missing(`origvar')
                    }
                }
                else {
                    qui gen double `xldate_`export_idx'' = `origvar' + 21916 if !missing(`origvar')
                }

                local new_export_varlist `new_export_varlist' `xldate_`export_idx''
                local date_converted_vars `date_converted_vars' `xldate_`export_idx''
            }
            else {
                local new_export_varlist `new_export_varlist' `var'
            }
        }
        if `has_date_cols' {
            local export_varlist `new_export_varlist'
        }
    }

    * Count variables and observations
    local nvars : word count `varlist'
    qui count `if' `in'
    local nobs = r(N)

    if `nvars' == 0 {
        di as error "cexport excel: no variables to export"
        exit 198
    }

    * Load plugin
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
            di as error "cexport excel: Could not load ctools plugin"
            exit 601
        }
    }

    * Build plugin arguments
    local opt_sheet = cond("`sheet'" != "", "sheet=`sheet'", "")
    local opt_replace = cond("`replace'" != "", "replace", "")
    local opt_nolabel = cond("`nolabel'" != "", "nolabel", "")
    local opt_verbose = cond("`verbose'" != "", "verbose", "")
    local opt_cell = cond(`"`cell'"' != "", `"cell=`cell'"', "")
    local opt_missing = cond(`"`missing'"' != "", `"missing=`missing'"', "")
    local opt_keepcellfmt = cond("`keepcellfmt'" != "", "keepcellfmt", "")

    * Pass variable names to plugin
    global CEXPORT_VARNAMES `varlist'

    * Pass variable types
    local vartypes ""
    foreach var of local export_varlist {
        local vtype : type `var'
        if "`vtype'" == "byte" {
            local vartypes `vartypes' 1
        }
        else if "`vtype'" == "int" {
            local vartypes `vartypes' 2
        }
        else if "`vtype'" == "long" {
            local vartypes `vartypes' 3
        }
        else if "`vtype'" == "float" {
            local vartypes `vartypes' 4
        }
        else if "`vtype'" == "double" {
            local vartypes `vartypes' 5
        }
        else {
            * String
            local vartypes `vartypes' 0
        }
    }
    global CEXPORT_VARTYPES `vartypes'

    * Build date column flags for plugin
    local date_cols ""
    local export_idx = 0
    foreach var of local export_varlist {
        local export_idx = `export_idx' + 1
        local origvar : word `export_idx' of `varlist'
        local fmt : format `origvar'
        local is_date = 0
        if substr("`fmt'", 1, 2) == "%t" | substr("`fmt'", 1, 3) == "%-t" {
            local is_date = 1
        }
        local date_cols `date_cols' `is_date'
    }
    global CEXPORT_DATE_COLS `date_cols'

    * Record start time
    timer clear 99
    timer on 99

    * Set string width metadata for flat buffer optimization
    _ctools_strw `export_varlist'

    * Call plugin with if `touse' so SF_ifobs() filters correctly
    capture noisily plugin call ctools_plugin `export_varlist' if `touse', ///
        "cexport_xlsx `using' `opt_sheet' `opt_firstrow' `opt_replace' `opt_nolabel' `opt_verbose' `opt_cell' `opt_missing' `opt_keepcellfmt'"

    local export_rc = _rc

    * Clean up
    macro drop CEXPORT_VARNAMES
    macro drop CEXPORT_VARTYPES
    macro drop CEXPORT_DATE_COLS

    if `export_rc' {
        di as error "Error exporting XLSX data (rc=`export_rc')"
        exit `export_rc'
    }

    timer off 99
    quietly timer list 99
    local elapsed = r(t99)

    * Display summary
    di as text ""
    di as text "file " as result `"`using'"' as text " saved"
    di as text "(" as result %12.0fc `nobs' as text " observations written)"

    if "`verbose'" != "" {
        di as text ""
        di as text "Time: " as result %9.3f `elapsed' as text " seconds"
    }

    * Return results
    return scalar N = `nobs'
    return scalar k = `nvars'
    return scalar time = `elapsed'
    return local filename `"`using'"'

    timer clear 99
end
