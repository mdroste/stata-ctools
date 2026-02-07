*! version 0.9.1 06Feb2026
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
    version 14.1

    * Parse the subcommand
    gettoken subcmd 0 : 0, parse(" ,")

    if "`subcmd'" == "excel" {
        * Dispatch to Excel import handler
        cimport_excel `0'
        return add
        exit
    }
    else if "`subcmd'" != "delimited" {
        di as error "cimport: unknown subcommand `subcmd'"
        di as error "Syntax: cimport delimited using filename [, options]"
        di as error "        cimport excel using filename [, options]"
        exit 198
    }

    * Now parse the rest with import delimited-style syntax
    * Support both "using filename" and just "filename" (like import delimited does)
    syntax [anything] [using/] [, Delimiters(string) VARNames(string) CLEAR ///
        CASE(string) ENCoding(string) BINDQuotes(string) ///
        STRIPQuotes ROWRange(string) COLRange(string) Verbose THReads(integer 0) ///
        ASFloat ASDOUBle NUMERICcols(numlist) STRINGcols(numlist) ///
        DECIMALSEParator(string) GROUPSEParator(string) ///
        MAXQUOTEDrows(integer 20) EMPTYlines(string) ///
        LOCale(string) PARSElocale]

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

    * Set default delimiter (auto-detect if not specified, matching Stata behavior)
    local plugin_delim "auto"
    if `"`delimiters'"' != "" {
        local plugin_delim `"`delimiters'"'
    }

    * Validate delimiter - for now only single character supported
    * For plugin: pass "tab" or "space" keyword instead of actual characters
    if `"`delimiters'"' != "" & length(`"`delimiters'"') != 1 {
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
    * Handle space delimiter - must pass as keyword since strtok uses space as separator
    if `"`delimiters'"' == " " {
        local plugin_delim "space"
    }

    * Parse varnames option
    * varnames(N) uses row N as header, skips rows 1 to N-1
    * varnames(nonames) means no header row
    local noheader = 0
    local headerrow = 1
    if "`varnames'" != "" {
        if "`varnames'" == "nonames" {
            local noheader = 1
            local headerrow = 0
        }
        else {
            * Must be a positive integer
            capture confirm integer number `varnames'
            if _rc != 0 | `varnames' < 1 {
                di as error "cimport: varnames() must be a positive integer or nonames"
                exit 198
            }
            local headerrow = `varnames'
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

    * Parse encoding - supports auto-detection and common encodings
    * Supported: utf-8, utf-16, iso-8859-1 (latin1), iso-8859-15 (latin9),
    *            windows-1252 (cp1252), ascii, macroman
    local encoding_opt ""
    if "`encoding'" != "" {
        local enc_lower = lower("`encoding'")
        * Normalize common encoding names
        if inlist("`enc_lower'", "utf-8", "utf8") {
            local encoding_opt "encoding=utf8"
        }
        else if inlist("`enc_lower'", "utf-16", "utf16", "utf-16le", "utf16le", "unicode") {
            local encoding_opt "encoding=utf16le"
        }
        else if inlist("`enc_lower'", "utf-16be", "utf16be") {
            local encoding_opt "encoding=utf16be"
        }
        else if inlist("`enc_lower'", "iso-8859-1", "iso88591", "latin1", "latin-1") {
            local encoding_opt "encoding=iso88591"
        }
        else if inlist("`enc_lower'", "iso-8859-15", "iso885915", "latin9", "latin-9") {
            local encoding_opt "encoding=iso885915"
        }
        else if inlist("`enc_lower'", "windows-1252", "windows1252", "cp1252", "win1252") {
            local encoding_opt "encoding=windows1252"
        }
        else if inlist("`enc_lower'", "ascii", "us-ascii") {
            local encoding_opt "encoding=ascii"
        }
        else if inlist("`enc_lower'", "macroman", "macintosh", "x-mac-roman") {
            local encoding_opt "encoding=macroman"
        }
        else {
            di as text "Warning: Unrecognized encoding '`encoding'', using auto-detection"
        }
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
    * NOTE: Stata's rowrange uses 1-based FILE line numbers (line 1 = header if present)
    * We convert to data row numbers for internal use
    local startrow = 0
    local endrow = 0
    local rowrange_header_only = 0
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

        * Convert from file line numbers to data row numbers
        * headerrow is the file line containing variable names (0 if nonames)
        * Data row 1 is at file line (headerrow + 1) if header exists, else line 1
        if `noheader' == 0 {
            * Has header at line `headerrow'
            * Check if rowrange only selects header/pre-header row(s)
            if `endrow' > 0 & `endrow' <= `headerrow' {
                local rowrange_header_only = 1
            }
            if `startrow' > 0 {
                local startrow = max(1, `startrow' - `headerrow')
            }
            if `endrow' > 0 {
                local endrow = `endrow' - `headerrow'
                if `endrow' < 1 local endrow = 0
            }
        }
    }

    * Parse colrange option
    local startcol = 0
    local endcol = 0
    if "`colrange'" != "" {
        * Parse start:end format
        local colonpos = strpos("`colrange'", ":")
        if `colonpos' > 0 {
            local startcol = substr("`colrange'", 1, `colonpos' - 1)
            local endcol = substr("`colrange'", `colonpos' + 1, .)
            if "`startcol'" == "" local startcol = 0
            if "`endcol'" == "" local endcol = 0
        }
        else {
            local startcol = `colrange'
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

    * Parse locale and parselocale options
    * parselocale enables locale-aware parsing of numbers
    * locale() specifies which locale to use; if not specified, uses system default
    * Common locale mappings:
    *   en_US, en_GB, en_AU, etc.: decimal=".", group=","
    *   de_DE, de_AT, nl_NL, etc.: decimal=",", group="."
    *   de_CH, fr_CH: decimal=".", group="'"
    *   fr_FR, it_IT, es_ES, pt_BR, etc.: decimal=",", group=" " or "."
    if "`parselocale'" != "" {
        * parselocale is specified - enable locale-aware parsing
        local locale_lower = lower("`locale'")

        * Map locale to decimal/group separators
        * If decimalseparator/groupseparator already specified, don't override
        if `"`decimalseparator'"' == "" {
            * Locales that use comma as decimal separator
            if inlist("`locale_lower'", "de_de", "de_at", "nl_nl", "nl_be", "fr_fr", "fr_be") | ///
               inlist("`locale_lower'", "it_it", "es_es", "pt_br", "pt_pt", "ru_ru", "pl_pl") | ///
               inlist("`locale_lower'", "cs_cz", "da_dk", "fi_fi", "el_gr", "hu_hu", "id_id") | ///
               inlist("`locale_lower'", "no_no", "ro_ro", "sk_sk", "sl_si", "sv_se", "tr_tr") | ///
               inlist("`locale_lower'", "uk_ua", "vi_vn", "bg_bg", "ca_es", "hr_hr", "et_ee") {
                local decimalsep ","
                * Most of these use period or space as group separator
                if `"`groupseparator'"' == "" {
                    * French and some others use space, most use period
                    if inlist("`locale_lower'", "fr_fr", "fr_be", "ru_ru", "pl_pl", "cs_cz") | ///
                       inlist("`locale_lower'", "fi_fi", "no_no", "sk_sk", "sv_se", "uk_ua") {
                        local groupsep " "
                    }
                    else {
                        local groupsep "."
                    }
                }
            }
            * Swiss locales use period decimal but apostrophe grouping
            else if inlist("`locale_lower'", "de_ch", "fr_ch", "it_ch") {
                local decimalsep "."
                if `"`groupseparator'"' == "" {
                    local groupsep "'"
                }
            }
            * English and similar locales use period decimal, comma grouping
            else if inlist("`locale_lower'", "en_us", "en_gb", "en_au", "en_ca", "en_nz") | ///
                    inlist("`locale_lower'", "en_ie", "en_za", "ja_jp", "ko_kr", "zh_cn") | ///
                    inlist("`locale_lower'", "zh_tw", "th_th", "ms_my", "en_in", "en_sg") | ///
                    "`locale_lower'" == "" {
                * Default: period decimal (already set above)
                if `"`groupseparator'"' == "" {
                    local groupsep ","
                }
            }
            * If locale not recognized, try to get system locale info
            else if "`locale_lower'" == "" {
                * No locale specified - try to detect system locale
                * Use period decimal as safe default
                local decimalsep "."
            }
            * Unrecognized locale - warn but continue with defaults
            else {
                di as text "Note: Unrecognized locale '`locale'', using default number format"
            }
        }
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
        di as text "Delimiter:  " as result cond(`"`delimiters'"' == "	", "tab", cond(`"`delimiters'"' == " ", "space", `"`delimiters'"'))
        di as text "Header row: " as result cond(`noheader' == 1, "no", "row `headerrow'")
        di as text "Case:       " as result "`case'"
        di as text "{hline 60}"
        di ""
    }

    * Build plugin arguments
    local opt_noheader = cond(`noheader' == 1, "noheader", "")
    local opt_headerrow = cond(`headerrow' > 1, "headerrow=`headerrow'", "")
    local opt_verbose = cond("`verbose'" != "", "verbose", "")
    local opt_case = "case=`case'"
    local opt_bindquotes = "bindquotes=`bindquotes'"

    * New options
    local opt_asfloat = cond("`asfloat'" != "", "asfloat", "")
    local opt_asdouble = cond("`asdouble'" != "", "asdouble", "")
    local opt_decimalsep = cond("`decimalsep'" != ".", "decimalsep=`decimalsep'", "")
    * Handle space groupsep specially (can't pass space character in space-delimited args)
    if "`groupsep'" == " " {
        local opt_groupsep "groupsep=space"
    }
    else {
        local opt_groupsep = cond("`groupsep'" != "", "groupsep=`groupsep'", "")
    }
    local opt_emptylines = cond("`emptylines'" != "skip", "emptylines=`emptylines'", "")
    local opt_maxquotedrows = cond(`maxquotedrows' != 20, "maxquotedrows=`maxquotedrows'", "")

    * Pass numericcols/stringcols via global macros (space-separated column numbers)
    if "`numcols_str'" != "" {
        global CIMPORT_NUMCOLS "`numcols_str'"
    }
    if "`strcols_str'" != "" {
        global CIMPORT_STRCOLS "`strcols_str'"
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
        "cimport `threads_code' scan `using' `plugin_delim' `opt_noheader' `opt_headerrow' `opt_verbose' `opt_bindquotes' `opt_asfloat' `opt_asdouble' `opt_decimalsep' `opt_groupsep' `opt_emptylines' `opt_maxquotedrows' `encoding_opt'"

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
    capture macro drop CIMPORT_NUMCOLS CIMPORT_STRCOLS

    if "`verbose'" != "" {
        di as text "  Found " as result `nobs' as text " rows, " as result `nvar' as text " columns"
    }

    * Handle empty file or file with header only - match Stata's behavior (rc=0, N=0, k=0)
    if `nvar' == 0 | `nobs' == 0 {
        timer off 11
        timer off 99
        quietly timer list 99
        local elapsed = r(t99)
        if `nobs' == 0 & `nvar' == 0 {
            di as text "(file is empty)"
        }
        else if `nobs' == 0 {
            di as text "(no data rows found)"
        }
        timer clear 11
        timer clear 99
        return scalar N = 0
        return scalar k = 0
        return scalar time = `elapsed'
        return local filename `"`using'"'
        exit 0
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

        * Try to create the variable; if name is reserved, fall back to v<i>
        local orig_vname "`vname'"
        local gen_ok = 1
        if `vtype' == 1 {
            * String variable
            if `vlen' < 1 local vlen = 1
            if `vlen' > 2045 local vlen = 2045
            capture quietly gen str`vlen' `vname' = ""
            if _rc != 0 {
                local gen_ok = 0
            }
        }
        else {
            * Numeric variable - use optimal storage type unless asfloat/asdouble specified
            if "`asfloat'" != "" {
                capture quietly gen float `vname' = .
            }
            else if "`asdouble'" != "" {
                capture quietly gen double `vname' = .
            }
            else if `ntype' == 4 {
                capture quietly gen byte `vname' = .
            }
            else if `ntype' == 3 {
                capture quietly gen int `vname' = .
            }
            else if `ntype' == 2 {
                capture quietly gen long `vname' = .
            }
            else if `ntype' == 1 {
                capture quietly gen float `vname' = .
            }
            else {
                capture quietly gen double `vname' = .
            }
            if _rc != 0 {
                local gen_ok = 0
            }
        }

        * If gen failed (reserved keyword, etc.), use fallback name v<i>
        if `gen_ok' == 0 {
            local vname "v`i'"
            * Ensure fallback name is unique
            capture confirm variable `vname'
            if !_rc {
                local suffix = 1
                local newname `vname'_`suffix'
                while (1) {
                    capture confirm variable `newname'
                    if _rc {
                        local vname `newname'
                        continue, break
                    }
                    local suffix = `suffix' + 1
                    local newname `vname'_`suffix'
                }
            }

            if `vtype' == 1 {
                quietly gen str`vlen' `vname' = ""
            }
            else {
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

            * Add variable label with original CSV header name
            label variable `vname' "`orig_vname'"
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
        "cimport `threads_code' load `using' `plugin_delim' `opt_noheader' `opt_headerrow' `opt_verbose' `opt_bindquotes' `opt_asfloat' `opt_asdouble' `opt_decimalsep' `opt_groupsep' `opt_emptylines' `opt_maxquotedrows' `encoding_opt'"

    local load_rc = _rc
    if `load_rc' {
        di as error "Error loading CSV data (rc=`load_rc')"
        exit `load_rc'
    }

    timer off 13

    * Apply rowrange filtering (post-import)
    * If rowrange only selected header row(s), drop all data
    if `rowrange_header_only' == 1 {
        quietly drop in 1/l
    }
    else if `startrow' > 0 | `endrow' > 0 {
        local first_keep = max(1, `startrow')
        if `endrow' > 0 {
            local last_keep = min(`endrow', _N)
        }
        else {
            local last_keep = _N
        }
        * Drop rows outside range (drop end first to preserve indices)
        if `last_keep' < _N {
            quietly drop in `=`last_keep'+1'/l
        }
        if `first_keep' > 1 {
            quietly drop in 1/`=`first_keep'-1'
        }
    }

    * Apply colrange filtering (post-import)
    if `startcol' > 0 | `endcol' > 0 {
        local first_col = max(1, `startcol')
        local total_cols = c(k)
        if `endcol' > 0 {
            local last_col = min(`endcol', `total_cols')
        }
        else {
            local last_col = `total_cols'
        }
        * Drop columns outside range
        if `last_col' < `total_cols' {
            forvalues i = `total_cols'(-1)`=`last_col'+1' {
                local vname : word `i' of `allvars'
                quietly drop `vname'
            }
        }
        if `first_col' > 1 {
            forvalues i = `=`first_col'-1'(-1)1 {
                local vname : word `i' of `allvars'
                quietly drop `vname'
            }
        }
    }

    timer off 99
    quietly timer list 99
    local elapsed = r(t99)

    * Display summary (matching import delimited output format)
    di as text "(" as result %12.0fc _N as text " observations read)"

    if "`verbose'" != "" {
        * Calculate Stata overhead
        capture local __plugin_time_total = _cimport_time_total
        if _rc != 0 local __plugin_time_total = 0
        local __stata_overhead = `elapsed' - `__plugin_time_total'
        if `__stata_overhead' < 0 local __stata_overhead = 0

        di as text ""
        di as text "{hline 55}"
        di as text "cimport timing breakdown:"
        di as text "{hline 55}"
        di as text "  C plugin internals:"
        di as text "    Memory map file:        " as result %8.4f _cimport_time_mmap " sec"
        di as text "    Parse CSV:              " as result %8.4f _cimport_time_parse " sec"
        di as text "    Type inference:         " as result %8.4f _cimport_time_infer " sec"
        di as text "    Cache conversion:       " as result %8.4f _cimport_time_cache " sec"
        di as text "    Store to Stata:         " as result %8.4f _cimport_time_store " sec"
        di as text "  {hline 53}"
        di as text "    C plugin total:         " as result %8.4f _cimport_time_total " sec"
        di as text "  {hline 53}"
        di as text "  Stata overhead:           " as result %8.4f `__stata_overhead' " sec"
        di as text "{hline 55}"
        di as text "    Wall clock total:       " as result %8.4f `elapsed' " sec"
        di as text "{hline 55}"

        * Throughput info
        if `elapsed' > 0 {
            tempname fh
            file open `fh' using `"`using'"', read binary
            file seek `fh' eof
            local fsize = r(loc)
            file close `fh'
            local mbps = (`fsize' / 1048576) / `elapsed'
            di as text "    Throughput:             " as result %9.1f `mbps' as text " MB/s"
        }

        * Display thread diagnostics
        capture local __threads_max = _cimport_threads_max
        if _rc == 0 {
            capture local __openmp_enabled = _cimport_openmp_enabled
            if _rc != 0 local __openmp_enabled = 0
            di as text ""
            di as text "  Thread diagnostics:"
            di as text "    OpenMP enabled:         " as result %8.0f `__openmp_enabled'
            di as text "    Max threads available:  " as result %8.0f `__threads_max'
            di as text "{hline 55}"
        }

        * Clean up timing scalars
        capture scalar drop _cimport_time_mmap _cimport_time_parse _cimport_time_infer
        capture scalar drop _cimport_time_cache _cimport_time_store _cimport_time_total
        capture scalar drop _cimport_threads_max _cimport_openmp_enabled
    }

    timer clear 11
    timer clear 12
    timer clear 13
    timer clear 99

    * Return results
    return scalar N = _N
    return scalar k = c(k)
    return scalar time = `elapsed'
    return local filename `"`using'"'

end

/*******************************************************************************
 * cimport_excel: Import Excel (.xlsx) files
 *
 * Syntax: cimport excel [using] filename [, options]
 *
 * Options:
 *   sheet(name)            - Worksheet to import (default: first sheet)
 *   cellrange([start][:end]) - Cell range to import (e.g., A1:D100, A1, :D100)
 *   firstrow               - First row contains variable names
 *   allstring              - Import all columns as strings
 *   case(preserve|lower|upper) - Variable name case handling
 *   clear                  - Clear current data before import
 *   verbose                - Display timing information
 ******************************************************************************/
program define cimport_excel, rclass
    version 14.1

    syntax [anything] [using/] [, SHEET(string) CELLRange(string) FIRSTrow ///
        ALLString CASE(string) CLEAR Verbose]

    * Handle filename - can come from using/ or as first positional argument
    if `"`using'"' == "" & `"`anything'"' != "" {
        local using `"`anything'"'
    }
    if `"`using'"' == "" {
        di as error "cimport excel: filename required"
        di as error "Syntax: cimport excel [using] filename [, options]"
        exit 198
    }

    * Validate file exists
    confirm file `"`using'"'

    * Validate .xlsx extension
    local ext = substr(`"`using'"', -5, 5)
    if lower("`ext'") != ".xlsx" {
        di as error "cimport excel: file must have .xlsx extension"
        exit 198
    }

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

    * Parse case option (default is preserve)
    if "`case'" != "" {
        if !inlist("`case'", "preserve", "lower", "upper") {
            di as error "cimport excel: case() must be preserve, lower, or upper"
            exit 198
        }
    }
    else {
        local case "preserve"
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
            di as error "cimport excel: Could not load ctools plugin"
            exit 601
        }
    }

    * Display info if verbose
    if "`verbose'" != "" {
        di as text ""
        di as text "{hline 60}"
        di as text "cimport excel: High-Performance Excel Import"
        di as text "{hline 60}"
        di as text "File:       " as result `"`using'"'
        if "`sheet'" != "" {
            di as text "Sheet:      " as result "`sheet'"
        }
        if "`cellrange'" != "" {
            di as text "Cell range: " as result "`cellrange'"
        }
        di as text "First row:  " as result cond("`firstrow'" != "", "variable names", "data")
        di as text "Case:       " as result "`case'"
        di as text "{hline 60}"
        di ""
    }

    * Build plugin arguments
    local opt_sheet = cond("`sheet'" != "", "sheet=`sheet'", "")
    local opt_cellrange = cond("`cellrange'" != "", "cellrange=`cellrange'", "")
    local opt_firstrow = cond("`firstrow'" != "", "firstrow", "")
    local opt_allstring = cond("`allstring'" != "", "allstring", "")
    local opt_verbose = cond("`verbose'" != "", "verbose", "")
    local opt_case = "case=`case'"

    * Record start time
    timer clear 99
    timer on 99

    * =========================================================================
    * PHASE 1: Scan the XLSX to get metadata
    * =========================================================================

    timer on 11

    if "`verbose'" != "" {
        di as text "Phase 1: Scanning XLSX file..."
    }

    capture noisily plugin call ctools_plugin, ///
        "cimport scan `using' `opt_sheet' `opt_cellrange' `opt_firstrow' `opt_allstring' `opt_case' `opt_verbose'"

    local scan_rc = _rc
    if `scan_rc' {
        di as error "Error scanning XLSX file (rc=`scan_rc')"
        exit `scan_rc'
    }

    * Retrieve metadata from global macros
    local nobs = ${_cimport_nobs}
    local nvar = ${_cimport_nvar}
    local varnames ${_cimport_varnames}
    local vartypes ${_cimport_vartypes}
    local numtypes ${_cimport_numtypes}
    local strlens ${_cimport_strlens}

    * Clean up global macros
    macro drop _cimport_nobs _cimport_nvar _cimport_varnames ///
               _cimport_vartypes _cimport_numtypes _cimport_strlens

    if "`verbose'" != "" {
        di as text "  Found " as result `nobs' as text " rows, " as result `nvar' as text " columns"
    }

    * Handle empty file
    if `nvar' == 0 | `nobs' == 0 {
        timer off 11
        timer off 99
        quietly timer list 99
        local elapsed = r(t99)
        if `nobs' == 0 & `nvar' == 0 {
            di as text "(file is empty)"
        }
        else if `nobs' == 0 {
            di as text "(no data rows found)"
        }
        timer clear 11
        timer clear 99
        return scalar N = 0
        return scalar k = 0
        return scalar time = `elapsed'
        return local filename `"`using'"'
        exit 0
    }
    timer off 11

    * =========================================================================
    * PHASE 2: Create variables in Stata
    * =========================================================================

    timer on 12
    if "`verbose'" != "" {
        di as text "Phase 2: Creating variables..."
    }

    quietly set obs 1

    * Create variables with appropriate types
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
            * Numeric variable
            if `ntype' == 4 {
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

    quietly set obs `nobs'

    timer off 12

    * =========================================================================
    * PHASE 3: Load data into variables
    * =========================================================================

    timer on 13

    if "`verbose'" != "" {
        di as text "Phase 3: Loading data..."
    }

    unab allvars : *

    capture noisily plugin call ctools_plugin `allvars', ///
        "cimport load `using' `opt_sheet' `opt_cellrange' `opt_firstrow' `opt_allstring' `opt_case' `opt_verbose'"

    local load_rc = _rc
    if `load_rc' {
        di as error "Error loading XLSX data (rc=`load_rc')"
        exit `load_rc'
    }

    timer off 13
    timer off 99
    quietly timer list 99
    local elapsed = r(t99)

    * Display summary
    di as text "(" as result %12.0fc _N as text " observations read)"

    if "`verbose'" != "" & `elapsed' > 0 {
        di as text ""
        di as text "Time:      " as result %9.3f `elapsed' as text " seconds"
    }

    timer clear 11
    timer clear 12
    timer clear 13
    timer clear 99

    * Return results
    return scalar N = _N
    return scalar k = c(k)
    return scalar time = `elapsed'
    return local filename `"`using'"'

end
