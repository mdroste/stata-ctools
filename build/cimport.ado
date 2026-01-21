*! version 1.1.0 21Jan2026
*! cimport: C-accelerated CSV import for Stata
*! Part of the ctools suite
*!
*! Description:
*!   High-performance replacement for import delimited using a C plugin
*!   with multi-threaded parallel parsing.
*!
*! Syntax:
*!   cimport delimited [using filename] [, options]
*!
*! Options:
*!   delimiters(string)  - Field delimiter(s) (default: ",")
*!   varnames(1|nonames) - First row is variable names (1) or data (nonames)
*!   clear               - Clear current data before import
*!   case(preserve|lower|upper) - Variable name case handling
*!   encoding(string)    - File encoding (currently only UTF-8 supported)
*!   bindquotes(strict|loose) - Quote handling mode
*!   stripquotes         - Remove surrounding quotes from string values
*!   rowrange([start][:end]) - Row range to import
*!   asfloat             - Import all numeric variables as float
*!   asdouble            - Import all numeric variables as double
*!   numericcols(numlist)  - Force specific columns to be numeric
*!   stringcols(numlist)   - Force specific columns to be string
*!   decimalseparator(char) - Decimal separator character (default: ".")
*!   groupseparator(char)  - Thousands grouping separator (default: none)
*!   maxquotedrows(#)    - Max rows to scan for quotes (default: 20)
*!   emptylines(skip|fill) - How to handle empty lines (default: skip)

program define cimport, rclass
    version 14.0

    * Parse the subcommand
    gettoken subcmd 0 : 0, parse(" ,")

    if "`subcmd'" != "delimited" {
        di as error "cimport: unknown subcommand `subcmd'"
        di as error "Syntax: cimport delimited using filename [, options]"
        exit 198
    }

    * Now parse the rest with import delimited-style syntax
    * Support both "using filename" and just "filename" (like import delimited does)
    syntax [anything] [using/] [, Delimiters(string) VARNames(string) CLEAR ///
        CASE(string) ENCoding(string) BINDQuotes(string) ///
        STRIPQuotes ROWRange(string) Verbose THReads(integer 0) ///
        ASFloat ASDOUBle NUMERICcols(numlist) STRINGcols(numlist) ///
        DECIMALSEParator(string) GROUPSEParator(string) ///
        MAXQUOTEDrows(integer 20) EMPTYlines(string)]

    * Handle filename - can come from using/ or as first positional argument
    if `"`using'"' == "" & `"`anything'"' != "" {
        * Filename provided without "using" keyword
        local using `"`anything'"'
    }
    if `"`using'"' == "" {
        di as error "cimport delimited: filename required"
        di as error "Syntax: cimport delimited [using] filename [, options]"
        exit 198
    }

    * Validate file exists
    confirm file `"`using'"'

    * Clear data if requested
    if "`clear'" != "" {
        clear
    }

    * Check that no data exists
    if _N > 0 | c(k) > 0 {
        di as error "data in memory would be lost"
        di as error "use the clear option to discard current data"
        exit 4
    }

    * Set default delimiter
    if `"`delimiters'"' == "" {
        local delimiters ","
    }

    * Validate delimiter - for now only single character supported
    * For plugin: pass "tab" keyword instead of actual tab character
    local plugin_delim `"`delimiters'"'
    if length(`"`delimiters'"') != 1 {
        * Handle special cases like tab
        if `"`delimiters'"' == "tab" | `"`delimiters'"' == "\t" {
            local delimiters "	"
            local plugin_delim "tab"
        }
        else {
            di as error "cimport: delimiter must be a single character"
            di as error "(multi-character delimiters not yet supported)"
            exit 198
        }
    }
    * Also check if delimiter is already a tab character
    if `"`delimiters'"' == "	" {
        local plugin_delim "tab"
    }

    * Parse varnames option
    local noheader = 0
    if "`varnames'" != "" {
        if "`varnames'" == "nonames" {
            local noheader = 1
        }
        else if "`varnames'" != "1" {
            di as error "cimport: varnames() must be 1 or nonames"
            exit 198
        }
    }

    * Parse case option (default is preserve)
    if "`case'" != "" {
        if !inlist("`case'", "preserve", "lower", "upper") {
            di as error "cimport: case() must be preserve, lower, or upper"
            exit 198
        }
    }
    else {
        local case "preserve"
    }

    * Parse encoding - currently only UTF-8 supported
    if "`encoding'" != "" & "`encoding'" != "utf-8" & "`encoding'" != "UTF-8" {
        di as text "Warning: cimport currently only supports UTF-8 encoding"
        di as text "         File will be read as UTF-8"
    }

    * Parse bindquotes - default is loose (matches Stata's import delimited default)
    if "`bindquotes'" != "" {
        if !inlist("`bindquotes'", "strict", "loose") {
            di as error "cimport: bindquotes() must be strict or loose"
            exit 198
        }
    }
    else {
        local bindquotes "loose"
    }

    * Parse rowrange option
    local startrow = 0
    local endrow = 0
    if "`rowrange'" != "" {
        * Parse start:end format
        local colonpos = strpos("`rowrange'", ":")
        if `colonpos' > 0 {
            local startrow = substr("`rowrange'", 1, `colonpos' - 1)
            local endrow = substr("`rowrange'", `colonpos' + 1, .)
            if "`startrow'" == "" local startrow = 0
            if "`endrow'" == "" local endrow = 0
        }
        else {
            local startrow = `rowrange'
        }
    }

    * Validate asfloat/asdouble - mutually exclusive
    if "`asfloat'" != "" & "`asdouble'" != "" {
        di as error "cimport: asfloat and asdouble are mutually exclusive"
        exit 198
    }

    * Parse decimalseparator - must be single character
    local decimalsep "."
    if `"`decimalseparator'"' != "" {
        if length(`"`decimalseparator'"') != 1 {
            di as error "cimport: decimalseparator() must be a single character"
            exit 198
        }
        local decimalsep `"`decimalseparator'"'
    }

    * Parse groupseparator - must be single character or empty
    local groupsep ""
    if `"`groupseparator'"' != "" {
        if length(`"`groupseparator'"') != 1 {
            di as error "cimport: groupseparator() must be a single character"
            exit 198
        }
        local groupsep `"`groupseparator'"'
    }

    * Parse emptylines option - default is skip
    if "`emptylines'" != "" {
        if !inlist("`emptylines'", "skip", "fill") {
            di as error "cimport: emptylines() must be skip or fill"
            exit 198
        }
    }
    else {
        local emptylines "skip"
    }

    * Validate maxquotedrows
    if `maxquotedrows' < 0 {
        di as error "cimport: maxquotedrows() must be non-negative"
        exit 198
    }

    * Convert numericcols/stringcols numlists to space-separated strings
    local numcols_str ""
    if "`numericcols'" != "" {
        local numcols_str "`numericcols'"
    }
    local strcols_str ""
    if "`stringcols'" != "" {
        local strcols_str "`stringcols'"
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
            di as error "cimport: Could not load ctools plugin"
            exit 601
        }
    }

    * Display info if verbose
    if "`verbose'" != "" {
        di as text ""
        di as text "{hline 60}"
        di as text "cimport delimited: High-Performance CSV Import"
        di as text "{hline 60}"
        di as text "File:       " as result `"`using'"'
        di as text "Delimiter:  " as result cond(`"`delimiters'"' == "	", "tab", `"`delimiters'"')
        di as text "Header row: " as result cond(`noheader' == 0, "yes (row 1)", "no")
        di as text "Case:       " as result "`case'"
        di as text "{hline 60}"
        di ""
    }

    * Build plugin arguments
    local opt_noheader = cond(`noheader' == 1, "noheader", "")
    local opt_verbose = cond("`verbose'" != "", "verbose", "")
    local opt_stripquotes = cond("`stripquotes'" != "", "stripquotes", "")
    local opt_case = "case=`case'"
    local opt_bindquotes = "bindquotes=`bindquotes'"

    * New options
    local opt_asfloat = cond("`asfloat'" != "", "asfloat", "")
    local opt_asdouble = cond("`asdouble'" != "", "asdouble", "")
    local opt_decimalsep = cond("`decimalsep'" != ".", "decimalsep=`decimalsep'", "")
    local opt_groupsep = cond("`groupsep'" != "", "groupsep=`groupsep'", "")
    local opt_emptylines = cond("`emptylines'" != "skip", "emptylines=`emptylines'", "")
    local opt_maxquotedrows = cond(`maxquotedrows' != 20, "maxquotedrows=`maxquotedrows'", "")

    * Pass numericcols/stringcols via global macros (space-separated column numbers)
    if "`numcols_str'" != "" {
        global _CIMPORT_NUMCOLS `numcols_str'
    }
    if "`strcols_str'" != "" {
        global _CIMPORT_STRCOLS `strcols_str'
    }

    * Build threads option string
    local threads_code ""
    if `threads' > 0 {
        local threads_code "threads(`threads')"
    }

    * Record start time
    timer clear 99
    timer on 99


    * =========================================================================
    * PHASE 1: Scan the CSV to get metadata
    * =========================================================================

    timer on 11

    if "`verbose'" != "" {
        di as text "Phase 1: Scanning CSV file..."
    }

    capture noisily plugin call ctools_plugin, ///
        "cimport `threads_code' scan `using' `plugin_delim' `opt_noheader' `opt_verbose' `opt_bindquotes' `opt_asfloat' `opt_asdouble' `opt_decimalsep' `opt_groupsep' `opt_emptylines' `opt_maxquotedrows'"

    local scan_rc = _rc
    if `scan_rc' {
        di as error "Error scanning CSV file (rc=`scan_rc')"
        exit `scan_rc'
    }

    * Retrieve metadata from global macros (SF_macro_save creates globals)
    local nobs = ${_cimport_nobs}
    local nvar = ${_cimport_nvar}
    local varnames ${_cimport_varnames}
    local vartypes ${_cimport_vartypes}
    local numtypes ${_cimport_numtypes}
    local strlens ${_cimport_strlens}

    * Clean up global macros
    macro drop _cimport_nobs _cimport_nvar _cimport_varnames ///
               _cimport_vartypes _cimport_numtypes _cimport_strlens
    capture macro drop _CIMPORT_NUMCOLS _CIMPORT_STRCOLS

    if "`verbose'" != "" {
        di as text "  Found " as result `nobs' as text " rows, " as result `nvar' as text " columns"
    }

    if `nobs' == 0 {
        di as error "No data rows found in CSV file"
        exit 2
    }
    timer off 11


    * =========================================================================
    * PHASE 2: Create variables in Stata
    * =========================================================================

    timer on 12
    if "`verbose'" != "" {
        di as text "Phase 2: Creating variables..."
    }

    * Set number of observations as 1 to create empty variables, change after
    quietly set obs 1

    * Create variables with appropriate types
    * Numeric subtypes: 0=double, 1=float, 2=long, 3=int, 4=byte
    local i = 1
    foreach vname of local varnames {
        local vtype : word `i' of `vartypes'
        local ntype : word `i' of `numtypes'
        local vlen : word `i' of `strlens'

        * Apply case transformation
        if "`case'" == "lower" {
            local vname = lower("`vname'")
        }
        else if "`case'" == "upper" {
            local vname = upper("`vname'")
        }

        * Make variable name unique if needed
        capture confirm variable `vname'
        if !_rc {
            local suffix = 1
            local newname `vname'`suffix'
            while (1) {
                capture confirm variable `newname'
                if _rc {
                    local vname `newname'
                    continue, break
                }
                local suffix = `suffix' + 1
                local newname `vname'`suffix'
            }
        }

        if `vtype' == 1 {
            * String variable
            if `vlen' < 1 local vlen = 1
            if `vlen' > 2045 local vlen = 2045
            quietly gen str`vlen' `vname' = ""
        }
        else {
            * Numeric variable - use optimal storage type unless asfloat/asdouble specified
            if "`asfloat'" != "" {
                quietly gen float `vname' = .
            }
            else if "`asdouble'" != "" {
                quietly gen double `vname' = .
            }
            else if `ntype' == 4 {
                quietly gen byte `vname' = .
            }
            else if `ntype' == 3 {
                quietly gen int `vname' = .
            }
            else if `ntype' == 2 {
                quietly gen long `vname' = .
            }
            else if `ntype' == 1 {
                quietly gen float `vname' = .
            }
            else {
                quietly gen double `vname' = .
            }
        }

        local i = `i' + 1
    }

    * Set number of observations
    quietly set obs `nobs'

    timer off 12

    * =========================================================================
    * PHASE 3: Load data into variables
    * =========================================================================

    timer on 13

    if "`verbose'" != "" {
        di as text "Phase 3: Loading data..."
    }

    * Get list of all variables we just created (in order)
    unab allvars : *

    capture noisily plugin call ctools_plugin `allvars', ///
        "cimport `threads_code' load `using' `plugin_delim' `opt_noheader' `opt_verbose' `opt_bindquotes' `opt_asfloat' `opt_asdouble' `opt_decimalsep' `opt_groupsep' `opt_emptylines' `opt_maxquotedrows'"

    local load_rc = _rc
    if `load_rc' {
        di as error "Error loading CSV data (rc=`load_rc')"
        exit `load_rc'
    }

    timer off 13

    timer off 99
    quietly timer list 99
    local elapsed = r(t99)

    * Display summary (matching import delimited output format)
    di as text "(" as result %12.0fc _N as text " observations read)"

    if "`verbose'" != "" & `elapsed' > 0 {
        di as text ""
        di as text "Time:      " as result %9.3f `elapsed' as text " seconds"
        * Get file size for throughput calculation
        tempname fh
        file open `fh' using `"`using'"', read binary
        file seek `fh' eof
        local fsize = r(loc)
        file close `fh'
        local mbps = (`fsize' / 1048576) / `elapsed'
        di as text "Speed:     " as result %9.1f `mbps' as text " MB/s"

        * Display thread diagnostics
        capture local __threads_max = _cimport_threads_max
        if _rc == 0 {
            capture local __openmp_enabled = _cimport_openmp_enabled
            if _rc != 0 local __openmp_enabled = 0
            di as text ""
            di as text "Thread diagnostics:"
            di as text "  OpenMP enabled:         " as result %4.0f `__openmp_enabled'
            di as text "  Max threads available:  " as result %4.0f `__threads_max'
        }
    }

    timer clear

    * Return results
    return scalar N = _N
    return scalar k = c(k)
    return scalar time = `elapsed'
    return local filename `"`using'"'

end
