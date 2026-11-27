*! version 1.3.2 20261124 github.com/mdroste/stata-ctools
*! original: Stata's export delimited (export_delimited.ado 1.0.9, June 13, 2020)
*! original: Stata's export excel (export_excel.ado 1.1.1, July 30, 2018)
*! original: see help cexport for details

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
    if inlist("`subcmd'", "delim", "delimi", "delimit", "delimite") local subcmd "delimited"
    if "`subcmd'" == "dbas" local subcmd "dbase"
    if inlist("`subcmd'", "exc", "exce") local subcmd "excel"

    if "`subcmd'" == "sasxport" & _caller()<16 local subcmd "sasxport5"
    if "`subcmd'" == "shp" {
        _cio_shp_export `0'
        return add
        exit
    }

    if "`subcmd'" == "dbase" {
        _cio_export `subcmd' `macval(0)'
        return add
        exit
    }
    if inlist("`subcmd'", "sas", "spss", "sasxport5", "sasxport8") {
        _cio_export `subcmd' `0'
        return add
        exit
    }

    if "`subcmd'" == "excel" {
        * Dispatch to Excel export handler
        cexport_excel `macval(0)'
        return add
        exit
    }
    else if "`subcmd'" != "delimited" {
        di as error "cexport: unknown subcommand `subcmd'"
        di as error "Supported file formats: delimited, excel, spss, sasxport5, sasxport8, dbase, shp, sas"
        exit 198
    }

    * Check that data is loaded (match export delimited rc=102)
    if c(k) == 0 {
        di as error "no variables defined"
        exit 102
    }

    * Parse the rest with export delimited-style syntax
    * Support both "using filename" and just "filename" (like export delimited does)
    * Abbreviations as export delimited takes them (q, novar)
    syntax [anything] [using/] [if] [in], [Delimiter(string) NOVARnames Quote ///
        NOQUOTEif REPLACE DATAfmt DATEString(string) NOLabel Verbose ///
        MMAP NOFSYNC DIRECT PREFAULT CRLF NOPARALLEL THReads(integer 0)]
    * The if is evaluated now, as export delimited does, before any helper
    * replaces the r() or e() results it may name.
    local iftouse ""
    if `"`if'"' != "" {
        tempvar iftouse
        mark `iftouse' `if' `in'
        local if "if `iftouse'"
    }

    * Handle filename - can come from using/ or as last positional argument
    if `"`using'"' == "" & `"`anything'"' != "" {
        gettoken using extra : anything
        if trim(`"`extra'"') != "" {
            di as error "specify using before the filename when selecting variables"
            exit 198
        }
        local anything ""
    }
    if `"`using'"' == "" {
        di as error "cexport delimited: filename required"
        di as error "Syntax: cexport delimited [varlist] [using] filename [, options]"
        exit 198
    }
    * ~ is the home directory, as for Stata's own file commands
    _ctools_home using `"`macval(using)'"'
    mata: st_local("__basename", pathbasename(st_local("using")))
    if !strpos(`"`__basename'"', ".") local using `"`using'.csv"'
    * Resolve varlist from anything, leaving out the sample marker.
    if strtrim(`"`anything'"') != "" {
        unab varlist : `anything'
    }
    else {
        qui ds
        local varlist `r(varlist)'
    }
    local varlist : list varlist - iftouse

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
        * A byte token survives whitespace and quote parsing in the dispatcher.
        mata: st_local("plugin_delim", "byte="+strofreal(ascii(st_local("delimiter"))))
    }

    * The writer displays value labels and the usual %t formats itself, read
    * straight from the numeric variable, where Stata would first build a
    * string copy of the whole column. fmtcheck asks the plugin whether it
    * displays a format; a plugin without it leaves everything to Stata.
    _ctools_load
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    local __cexport_fmtq "%td"
    local __cexport_fmtok 0
    capture plugin call ctools_plugin, "cexport fmtcheck"
    local __cx_cells = (_rc == 0 & "`__cexport_fmtok'" == "1")
    capture confirm number 0

    * Handle value labels: decode labeled variables unless nolabel specified
    local export_varlist ""
    local decoded_vars ""
    local original_names ""
    * Preformatted string columns, which the writer never quotes
    local formatted_vars ""
    local ndec 0
    if "`nolabel'" == "" {
        * Default: export value labels
        local __cx_i = 0
        foreach var of local varlist {
            local ++__cx_i
            * Check if variable has a value label attached
            local vallbl : value label `var'
            * The writer looks the label text up itself when the label fits
            * its mapping; a labeled date, a label that does not fit and a
            * label never defined are decoded below.
            local __cx_ok 0
            if "`vallbl'" != "" & `__cx_cells' {
                local labfmt : format `var'
                if !(substr("`labfmt'",1,2) == "%t" | substr("`labfmt'",1,3) == "%-t") {
                    mata: _cexport_pack_label(st_local("vallbl"), `__cx_i')
                }
            }
            if `__cx_ok' {
                local __cexport_cell_`__cx_i' 2
                local export_varlist `export_varlist' `var'
            }
            else if "`vallbl'" != "" {
                * Has value label - decode to temporary string variable. The
                * tempvar is named by a counter, so any variable name fits; the
                * replace statements below widen the decoded str# as needed.
                local ++ndec
                tempvar dec`ndec'
                local dname `dec`ndec''
                * A label attached but never defined labels nothing: as
                * natively, every value is then written as a number.
                mata: st_local("__lblexists", strofreal(st_vlexists(st_local("vallbl"))))
                if `__lblexists' qui decode `var', gen(`dname')
                else qui gen str1 `dname' = ""
                local labfmt : format `var'
                * Unlabeled values are written as in an unlabeled column: 16
                * significant digits for doubles (%18.0g) and 8 for floats
                * (%14.0g of the value rounded to float precision), whatever
                * the display format; dates in their own format.
                local labfloat 0
                if !(substr("`labfmt'",1,2) == "%t" | substr("`labfmt'",1,3) == "%-t") {
                    local labfmt "%18.0g"
                    if "`: type `var''" == "float" local labfloat 1
                }
                * Under quote, native export quotes only the label text; a value
                * without a label is still written as an unquoted number (and
                * system missing as an empty unquoted field). Quote the label
                * cells here and mark the column preformatted so the writer
                * leaves the numeric fallbacks bare.
                if "`quote'" != "" {
                    qui replace `dname' = char(34) + subinstr(`dname', char(34), char(34)+char(34), .) + char(34) if `dname'!=""
                    local formatted_vars `formatted_vars' `dname'
                }
                if `labfloat' qui replace `dname' = strtrim(string(real(string(`var',"%15.7e")),"%14.0g")) if `dname'=="" & `var' != .
                else qui replace `dname' = strtrim(string(`var',"`labfmt'")) if `dname'=="" & `var' != .
                local export_varlist `export_varlist' `dname'
                local decoded_vars `decoded_vars' `dname'
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

    * Prepare display-format fields; the C writer still owns CSV encoding and I/O.
    * Native export writes date/time variables (%t and %-t) with their display
    * format by default; datafmt extends that to every numeric variable.
    local new_export_varlist ""
    local export_idx = 0
    foreach var of local export_varlist {
        local ++export_idx
        local origvar : word `export_idx' of `varlist'
        local fmt : format `origvar'
        local is_time = substr("`fmt'",1,2) == "%t" | substr("`fmt'",1,3) == "%-t"
        capture confirm numeric variable `var'
        if !_rc & "`__cexport_cell_`export_idx''" == "2" {
            * its value label is displayed by the writer
            local new_export_varlist `new_export_varlist' `var'
        }
        else if !_rc & ("`datafmt'" != "" | `is_time') {
            if `"`datestring'"' != "" & `is_time' local fmt `"`datestring'"'
            * A %t format the writer displays, on values Stata displays as
            * dates (years 100 to 9999; others it shows as numbers), needs no
            * formatted copy of the column.
            local __cx_ok 0
            if `is_time' & `__cx_cells' {
                local __cexport_fmtq `"`fmt'"'
                local __cexport_fmtok 0
                capture plugin call ctools_plugin, "cexport fmtcheck"
                if !_rc & "`__cexport_fmtok'" == "1" {
                    local __cx_unit = substr(`"`fmt'"', strpos(`"`fmt'"', "t") + 1, 1)
                    local __cx_lo .
                    if "`__cx_unit'" == "d" {
                        local __cx_lo -679350
                        local __cx_hi 2936549
                    }
                    else if "`__cx_unit'" == "c" {
                        local __cx_lo -58695840000000
                        local __cx_hi 253717919999999
                    }
                    else if "`__cx_unit'" == "w" {
                        local __cx_lo -96720
                        local __cx_hi 418079
                    }
                    else if "`__cx_unit'" == "m" {
                        local __cx_lo -22320
                        local __cx_hi 96479
                    }
                    else if "`__cx_unit'" == "q" {
                        local __cx_lo -7440
                        local __cx_hi 32159
                    }
                    else if "`__cx_unit'" == "h" {
                        local __cx_lo -3720
                        local __cx_hi 16079
                    }
                    else if "`__cx_unit'" == "y" {
                        local __cx_lo 100
                        local __cx_hi 9999
                    }
                    if `__cx_lo' < . {
                        quietly summarize `var', meanonly
                        if r(N) == 0 local __cx_ok 1
                        else if floor(r(min)) >= `__cx_lo' & floor(r(max)) <= `__cx_hi' local __cx_ok 1
                    }
                }
                capture confirm number 0
            }
            if `__cx_ok' {
                local __cexport_cell_`export_idx' 1
                local __cexport_fmt_`export_idx' `"`fmt'"'
                local __cexport_fmt_`export_idx'_len = strlen(`"`fmt'"')
                local new_export_varlist `new_export_varlist' `var'
                continue
            }
            * As export delimited does, a decimal-comma format is written with
            * a decimal point (%9,2fc as %9.2fc) and only leading blanks
            * removed; any other comma format loses its thousands separators
            * and zero padding (%012.0gc as %12.0g), and a left-justified one
            * is not zero padded (%-08.2f).
            local trim strtrim
            if !`is_time' {
                if strpos("`fmt'", ",") {
                    local fmt = subinstr("`fmt'", ",", ".", 1)
                    local trim strltrim
                }
                else {
                    if substr("`fmt'", 1, 3) == "%-0" local fmt = "%-" + substr("`fmt'", 4, .)
                    if substr("`fmt'", -1, 1) == "c" {
                        local fmt = substr("`fmt'", 1, length("`fmt'") - 1)
                        if substr("`fmt'", 1, 2) == "%0" local fmt = "%" + substr("`fmt'", 3, .)
                    }
                }
            }
            tempvar formatted
            quietly gen `formatted' = cond(`var' == ., "", `trim'(string(`var', `"`fmt'"')))
            local new_export_varlist `new_export_varlist' `formatted'
            local formatted_vars `formatted_vars' `formatted'
        }
        else local new_export_varlist `new_export_varlist' `var'
    }
    local export_varlist `new_export_varlist'

    * Count variables and observations
    local nvars : word count `varlist'
    qui count `if' `in'
    local nobs = r(N)

    if `nvars' == 0 {
        di as error "cexport: no variables to export"
        exit 198
    }

    * Load the platform-appropriate ctools plugin if not already loaded
    _ctools_load
    * Stata scopes plugin registrations to the calling ado program.
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    capture confirm number 0

    * Build plugin arguments
    local opt_noheader = cond("`novarnames'" != "", "noheader", "")
    local opt_quote = cond("`quote'" != "", "quote", "")
    local opt_noquoteif = cond("`noquoteif'" != "", "noquoteif", "")
    local opt_verbose = cond("`verbose'" != "", "verbose", "")

    * Performance options. nofsync is accepted for compatibility and has no
    * effect: cexport no longer waits for the file to reach the disk.
    local opt_mmap = cond("`mmap'" != "", "mmap", "")
    local opt_direct = cond("`direct'" != "", "direct", "")
    local opt_prefault = cond("`prefault'" != "", "prefault", "")
    local opt_crlf = cond("`crlf'" != "", "crlf", "")
    local opt_noparallel = cond("`noparallel'" != "", "noparallel", "")

    * Build threads option string
    local threads_code ""
    if `threads' > 0 {
        local threads_code "threads(`threads')"
    }



    * Pass variable storage types to the plugin
    * Types: 1=byte, 2=int, 3=long, 4=float, 5=double, 0=string
    local vartypes ""
    foreach var of local export_varlist {
        local vtype : type `var'
        if `: list var in formatted_vars' {
            local vartypes `vartypes' 6
        }
        else if "`vtype'" == "byte" {
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


    * Each metadata field has an independent local; no aggregate buffer limit.
    local metadata_count : word count `varlist'
    forvalues metadata_i=1/`metadata_count' {
        local __cexport_name_`metadata_i' : word `metadata_i' of `varlist'
        if "`firstrow'" == "varlabels" {
            local header_var : word `metadata_i' of `varlist'
            * Labels may come from someone else's file: Mata only
            mata: st_local("__cexport_name_`metadata_i'", st_varlabel("`header_var'") != "" ? st_varlabel("`header_var'") : st_local("__cexport_name_`metadata_i'"))
        }
        local __cexport_type_`metadata_i' : word `metadata_i' of `vartypes'
        local __cexport_date_`metadata_i' : word `metadata_i' of `date_cols'
        if "`__cexport_date_`metadata_i''" == "" local __cexport_date_`metadata_i' 0
    }
    local __cexport_filename `"`using'"'
    local __cexport_filename_len = strlen(`"`__cexport_filename'"')
    local __cexport_sheet `"`sheet'"'
    local __cexport_sheet_len = strlen(`"`__cexport_sheet'"')
    local __cexport_missing `"`missing'"'
    local __cexport_missing_len = strlen(`"`__cexport_missing'"')
    local opt_replace = cond("`replace'" != "", "replace", "")

    * Record start time
    local __ct_elapsed99 = 0
    if "`__ct_elapsed99'" == "" local __ct_elapsed99 = 0
    quietly _ctools_clock __ct_start99

    * Set string width metadata for flat buffer optimization
    _ctools_strw `export_varlist'

    * Stata answers an if one observation at a time. Without one the plugin
    * skips the check; the if marked after syntax is passed last, as a
    * sample marker the plugin reads in parallel.
    local opt_select "skipif"
    local marker ""
    local ifcond ""
    if `"`if'"' != "" {
        local marker "`iftouse'"
        local opt_select "touse"
    }

    * Call the C plugin
    * Plugin expects options only; paths/metadata are length-checked locals.
    * Use export_varlist (may contain decoded temp vars for value labels)
    capture noisily plugin call ctools_plugin `export_varlist' `marker' `ifcond' `in', ///
        "cexport `threads_code' `plugin_delim' `opt_replace' `opt_noheader' `opt_quote' `opt_noquoteif' `opt_verbose' `opt_mmap' `opt_direct' `opt_prefault' `opt_crlf' `opt_noparallel' `opt_select'"

    local export_rc = _rc

    * Clean up global macros

    * The plugin reports a write-protected destination as r(608)
    if `export_rc' == 608 {
        di as error `"file `using' is read-only; cannot be modified or erased"'
        exit 603
    }
    if `export_rc' {
        di as error "Error exporting data (rc=`export_rc')"
        exit `export_rc'
    }

    quietly _ctools_clock __ct_now

    local __ct_elapsed99 = `__ct_elapsed99' + (`__ct_now' - `__ct_start99')
    * Elapsed seconds are held in caller-local timing state.
    local elapsed = `__ct_elapsed99'

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

    local __ct_elapsed99 = 0
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
    * Abbreviations as export excel takes them (sh(), date(), miss())
    syntax [anything] [using/] [if] [in], [SHeet(string asis) FIRSTrow(string) ///
        REPLACE DATAfmt DATEstring(string) NOLabel Verbose ///
        CELL(string) MISSing(string) KEEPCELLfmt THReads(integer 0) ///
        SHEETMODify SHEETREPlace LOCALE(string)]
    * As export excel does, missing() holds at most 16 bytes
    if strlen(`"`macval(missing)'"') > 16 {
        di as error "invalid string in option {bf:missing()}"
        di as error "{p 4 4 2}The string is too long.{p_end}"
        exit 198
    }
    * The if is evaluated now, as export excel does, before any helper
    * replaces the r() or e() results it may name.
    local iftouse ""
    if `"`if'"' != "" {
        tempvar iftouse
        mark `iftouse' `if' `in'
        local if "if `iftouse'"
    }

    _cexport_excel_sheet `macval(sheet)'
    mata: st_local("sheet", st_global("r(sheet)"))
    local sheet_action `r(action)'
    * As export excel does: the separate sheetmodify and sheetreplace options
    * exclude each other and take precedence over the sheet() suboption.
    local sheet_option `"sheet(..., `=substr("`sheet_action'", 6, .)')"'
    if "`sheetmodify'" != "" & "`sheetreplace'" != "" {
        di as error "invalid syntax; option sheetmodify may not be combined with option sheetreplace"
        exit 198
    }
    if "`sheetmodify'`sheetreplace'" != "" {
        local sheet_action "`sheetmodify'`sheetreplace'"
        local sheet_option "`sheet_action'"
    }
    if "`sheet_action'" != "" & "`replace'" != "" {
        di as error "invalid syntax; option `sheet_option' may not be combined with option replace"
        exit 198
    }

    * Handle filename - can come from using/ or as last positional argument
    if `"`macval(using)'"' == "" & `"`macval(anything)'"' != "" {
        local __tokens `"`macval(anything)'"'
        local __filename ""
        local __varlist ""
        while `"`macval(__tokens)'"' != "" {
            gettoken __tok __tokens : __tokens
            * A one-line if reparses its command and expands the text again.
            if `"`macval(__tokens)'"' == "" {
                mata: st_local("__filename", st_local("__tok"))
            }
            else {
                local __varlist `"`macval(__varlist)' `macval(__tok)'"'
            }
        }
        if `"`macval(__filename)'"' == "" {
            di as error "cexport excel: filename required"
            di as error "Syntax: cexport excel [varlist] [using] filename [, options]"
            exit 198
        }
        local using `"`macval(__filename)'"'
        local anything `"`macval(__varlist)'"'
    }
    if `"`macval(using)'"' == "" {
        di as error "cexport excel: filename required"
        di as error "Syntax: cexport excel [varlist] [using] filename [, options]"
        exit 198
    }
    * Resolve varlist from anything.
    * Must happen BEFORE marksample, which creates a touse tempvar
    * that ds would include in the variable list.
    if strtrim(`"`macval(anything)'"') != "" {
        unab varlist : `macval(anything)'
    }
    else {
        qui ds
        local varlist `r(varlist)'
    }
    local varlist : list varlist - iftouse

    * Mark sample
    marksample touse, novarlist

    * Validate cell option format (e.g., A1, B5, AA10)
    if `"`macval(cell)'"' != "" {
        local cell_upper = upper(`"`macval(cell)'"')
        * Check format: 1-3 letters followed by 1 or more digits
        if !regexm("`cell_upper'", "^[A-Z]+[0-9]+$") {
            di as error "cexport excel: cell() must be a valid cell reference (e.g., A1, B5)"
            exit 198
        }
        local cell `"`cell_upper'"'
    }

    * ~ is the home directory, as for Stata's own file commands
    _ctools_home using `"`macval(using)'"'
    * As export excel does: a file name without an extension gets .xlsx; one
    * ending in .xls is a binary workbook; any other name (.xlsm, .xltx, .dat)
    * is kept as given and holds an Office Open XML workbook.
    mata: st_local("__xlbase", pathbasename(st_local("using")))
    if !strpos(`"`macval(__xlbase)'"', ".") {
        local using `"`macval(using)'.xlsx"'
    }
    local excel_engine "cexport_xlsx"
    if lower(substr(`"`macval(using)'"',-4,4)) == ".xls" local excel_engine "cexport_xls"

    * Parse firstrow option
    local opt_firstrow = "nofirstrow"
    if "`macval(firstrow)'" != "" & substr("variables",1,strlen("`macval(firstrow)'")) == "`macval(firstrow)'" local firstrow "variables"
    if strlen("`macval(firstrow)'") >= 4 & substr("varlabels",1,strlen("`macval(firstrow)'")) == "`macval(firstrow)'" local firstrow "varlabels"
    if "`macval(firstrow)'" != "" {
        if "`macval(firstrow)'" == "nonames" {
            local opt_firstrow = "nofirstrow"
        }
        else if inlist("`macval(firstrow)'", "variables", "varlabels") {
            local opt_firstrow "firstrow"
        }
        else {
            di as error "cexport excel: firstrow() must be 'variables' or 'nonames'"
            exit 198
        }
    }

    * Handle value labels
    local export_varlist ""
    local decoded_vars ""
    local ndec 0
    if "`nolabel'" == "" {
        foreach var of local varlist {
            local vallbl : value label `var'
            if "`vallbl'" != "" {
                * A counter names the tempvar, so any variable name fits
                local ++ndec
                tempvar dec`ndec'
                local dname `dec`ndec''
                * A label attached but never defined labels nothing (native
                * then writes every value as below).
                mata: st_local("__lblexists", strofreal(st_vlexists(st_local("vallbl"))))
                if `__lblexists' qui decode `var', gen(`dname')
                else qui gen str1 `dname' = ""
                * Unlabeled values, extended missing values included (.b), are
                * written as text in %9.0g, whatever the storage type or
                * display format, as native export excel does; only system
                * missing is left empty
                qui replace `dname' = string(`var',"%9.0g") if `dname'=="" & `var' != .
                local export_varlist `export_varlist' `dname'
                local decoded_vars `decoded_vars' `dname'
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
    if "`datafmt'" != "" | `"`macval(datestring)'"' != "" {
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
            * The old %d is a daily format, as export excel reads it
            if substr("`fmt'", 1, 2) == "%d" | substr("`fmt'", 1, 3) == "%-d" local is_date_fmt = 1

            if `is_date_fmt' {
                * This is a date/time variable - convert to formatted string
                tempvar fmtvar_`export_idx'

                * Use datestring format if specified, otherwise use variable's display format
                * Note: missing values should become empty strings (like Stata's behavior)
                if `"`macval(datestring)'"' != "" {
                    qui gen str `fmtvar_`export_idx'' = cond(missing(`origvar'), "", string(`origvar', `"`macval(datestring)'"'))
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
            * The old %d is a daily format, as export excel reads it
            if substr("`fmt'", 1, 2) == "%d" | substr("`fmt'", 1, 3) == "%-d" local is_date_fmt = 1
            if substr("`fmt'", 1, 2) == "%d" local fmt = "%t" + substr("`fmt'", 2, .)
            if substr("`fmt'", 1, 3) == "%-d" local fmt = "%-t" + substr("`fmt'", 3, .)
            * Determine date type: strip % and optional -
            local fmtsub = substr("`fmt'", 2, .)
            if substr("`fmtsub'", 1, 1) == "-" local fmtsub = substr("`fmtsub'", 2, .)
            local datetype = substr("`fmtsub'", 2, 1)
            * A generic %tg value has no calendar unit: write the raw number
            if `is_date_fmt' & !inlist("`datetype'", "d", "m", "q", "w", "h", "y", "b", "c", "C") local is_date_fmt = 0

            if `is_date_fmt' {
                local has_date_cols = 1

                tempvar xldate_`export_idx'
                if "`datetype'" == "b" {
                    * Business dates count calendar days: %tb<cal>[:display]
                    local bcal = substr("`fmtsub'", 3, .)
                    if strpos("`bcal'", ":") local bcal = substr("`bcal'", 1, strpos("`bcal'", ":") - 1)
                    qui gen double `xldate_`export_idx'' = dofb(`origvar', "`bcal'") + 21916 if !missing(`origvar')
                }
                else if "`datetype'" == "d" {
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
                        qui gen double `xldate_`export_idx'' = cofC(`origvar') / 86400000 + 21916 if !missing(`origvar')
                    }
                    else {
                        qui gen double `xldate_`export_idx'' = `origvar' / 86400000 + 21916 if !missing(`origvar')
                    }
                }
                else {
                    qui gen double `xldate_`export_idx'' = `origvar' + 21916 if !missing(`origvar')
                }

                * Excel inserts a nonexistent 29feb1900 at serial 60.
                qui replace `xldate_`export_idx'' = `xldate_`export_idx'' - 1 if `xldate_`export_idx'' < 61
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
    * As export excel does: an empty dataset is r(2000), an empty selection
    * r(198)
    if _N == 0 {
        di as error "no observations"
        exit 2000
    }
    if `nobs' == 0 {
        if "`excel_engine'" == "cexport_xls" di as error "too many or no observations specified"
        else di as error "observations must be between 1 and 1048576"
        exit 198
    }

    * Load plugin
    _ctools_load
    * Stata scopes plugin registrations to the calling ado program.
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    capture confirm number 0

    * Build plugin arguments
    local opt_replace = cond("`replace'" != "", "replace", "")
    local opt_nolabel = cond("`nolabel'" != "", "nolabel", "")
    local opt_verbose = cond("`verbose'" != "", "verbose", "")
    local opt_cell = cond(`"`macval(cell)'"' != "", `"cell=`macval(cell)'"', "")
    local opt_keepcellfmt = cond("`keepcellfmt'" != "", "keepcellfmt", "")

    * Pass variable names to plugin


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


    * Cell style codes for the plugin, assigned from the display format as
    * export excel assigns them: 1 date (m/d/yy), 3 date and time
    * (m/d/yy h:mm), 4 time of day (h:mm:ss AM/PM) when a %tc format starts
    * with an hour, minute or am/pm code, 2 a %tc format that shows seconds
    * with s (native writes empty text; cexport writes the value with
    * milliseconds), 5 0 for %tg; numbers 7 0.00 (.2f), 9 #,##0.00 (.2fc),
    * 8 #,##0 (gc), 10 0.00E+00 (.2e) and 6 0 otherwise. 0 is General: text,
    * and numbers under keepcellfmt, which keeps a cell's format for
    * everything but %t formats.
    local date_cols ""
    local export_idx = 0
    foreach var of local export_varlist {
        local export_idx = `export_idx' + 1
        local origvar : word `export_idx' of `varlist'
        local fmt : format `origvar'
        if substr("`fmt'", 1, 2) == "%d" local fmt = "%t" + substr("`fmt'", 2, .)
        if substr("`fmt'", 1, 3) == "%-d" local fmt = "%-t" + substr("`fmt'", 3, .)
        local is_date = 0
        capture confirm numeric variable `var'
        local isnum = !_rc
        if `isnum' & (substr("`fmt'", 1, 2) == "%t" | substr("`fmt'", 1, 3) == "%-t") {
            * The unit follows %t or %-t; %tg is written as its raw number
            local unit = substr(subinstr("`fmt'", "%-t", "%t", 1), 3, 1)
            local custom = substr(subinstr("`fmt'", "%-t", "%t", 1), 4, .)
            local is_date = 5
            if inlist("`unit'", "d", "m", "q", "w", "h", "y", "b") local is_date = 1
            else if inlist("`unit'", "c", "C") {
                local is_date = 3
                if inlist(substr("`custom'", 1, 1), "H", "h", "a") | inlist(substr("`custom'", 1, 2), "MM", "mm") local is_date = 4
                if strpos("`custom'", "s") local is_date = 2
            }
        }
        if `isnum' & !`is_date' & "`keepcellfmt'" == "" {
            local is_date = 6
            if substr("`fmt'", -2, 2) == "gc" local is_date = 8
            else if substr("`fmt'", -4, 4) == ".2fc" local is_date = 9
            else if substr("`fmt'", -3, 3) == ".2f" local is_date = 7
            else if substr("`fmt'", -3, 3) == ".2e" local is_date = 10
        }
        local date_cols `date_cols' `is_date'
    }


    * Each metadata field has an independent local; no aggregate buffer limit.
    local metadata_count : word count `varlist'
    forvalues metadata_i=1/`metadata_count' {
        local __cexport_name_`metadata_i' : word `metadata_i' of `varlist'
        if "`macval(firstrow)'" == "varlabels" {
            local header_var : word `metadata_i' of `varlist'
            * Labels may come from someone else's file: Mata only
            mata: st_local("__cexport_name_`metadata_i'", st_varlabel("`header_var'") != "" ? st_varlabel("`header_var'") : st_local("__cexport_name_`metadata_i'"))
        }
        local __cexport_type_`metadata_i' : word `metadata_i' of `vartypes'
        local __cexport_date_`metadata_i' : word `metadata_i' of `date_cols'
        if "`__cexport_date_`metadata_i''" == "" local __cexport_date_`metadata_i' 0
    }
    * Pass literal text directly to the plugin, without macro re-expansion.
    foreach field in filename sheet missing {
        local source = cond("`field'" == "filename", "using", "`field'")
        mata: st_local("__cexport_`field'", st_local("`source'"))
        mata: st_local("__cexport_`field'_len", strofreal(strlen(st_local("`source'"))))
    }
    local opt_replace = cond("`replace'" != "", "replace", "")

    * Record start time
    local __ct_elapsed99 = 0
    if "`__ct_elapsed99'" == "" local __ct_elapsed99 = 0
    quietly _ctools_clock __ct_start99

    * Set string width metadata for flat buffer optimization
    _ctools_strw `export_varlist'

    local threads_code ""
    if `threads' > 0 local threads_code "threads(`threads')"

    * Call plugin with if `touse' so SF_ifobs() filters correctly
    capture noisily plugin call ctools_plugin `export_varlist' if `touse', ///
        "`excel_engine' `threads_code' `opt_firstrow' `opt_replace' `opt_nolabel' `opt_verbose' `opt_cell' `opt_keepcellfmt' `sheet_action'"

    local export_rc = _rc

    * Clean up

    * The plugin reports a write-protected destination as r(608)
    if `export_rc' == 608 {
        di as error `"file `macval(using)' could not be saved"'
        exit 603
    }
    if `export_rc' {
        di as error "Error exporting XLSX data (rc=`export_rc')"
        exit `export_rc'
    }

    quietly _ctools_clock __ct_now

    local __ct_elapsed99 = `__ct_elapsed99' + (`__ct_now' - `__ct_start99')
    * Elapsed seconds are held in caller-local timing state.
    local elapsed = `__ct_elapsed99'

    * Display summary
    di as text ""
    di as text "file " as result `"`macval(using)'"' as text " saved"
    di as text "(" as result %12.0fc `nobs' as text " observations written)"

    if "`verbose'" != "" {
        di as text ""
        di as text "Time: " as result %9.3f `elapsed' as text " seconds"
    }

    * Return results
    return scalar N = `nobs'
    return scalar k = `nvars'
    return scalar time = `elapsed'
    return local filename `"`macval(using)'"'

    local __ct_elapsed99 = 0
end

* Parse the nested native sheet(name, modify|replace) syntax independently.
* As export excel does: one name, quoted or not, then modify or replace
* spelled in full; when both are given the first applies.
program define _cexport_excel_sheet, rclass
    version 14.1
    gettoken sheet rest : 0, parse(" ,") quotes
    if `"`macval(sheet)'"' == "," {
        local rest `", `macval(rest)'"'
        local sheet ""
    }
    else gettoken sheet : sheet
    gettoken comma rest : rest, parse(" ,")
    local action ""
    local bad = !inlist(`"`comma'"', "", ",") | (`"`comma'"' == "" & `"`macval(rest)'"' != "")
    while !`bad' & `"`macval(rest)'"' != "" {
        gettoken word rest : rest, quotes
        if !inlist(`"`macval(word)'"', "modify", "replace") local bad 1
        else if "`action'" == "" local action "sheet`word'"
    }
    if `bad' {
        di as error "invalid argument in option sheet()"
        exit 198
    }
    return local sheet `"`macval(sheet)'"'
    return local action "`action'"
end

version 14.1
mata:
/* The value label named label in cdecode's mapping format, for the writer:
   "<count> " then "<value> <bytes> <text> " per mapping and a closing "|".
   Sets the caller's locals __cexport_map_<j>, its _len and
   __cexport_mapwidth_<j> (the longest text) and __cx_ok to 1; leaves __cx_ok
   0 for a label that is not defined, has no mappings, holds a text over
   2,045 bytes or a NUL byte, or does not fit a macro. */
void _cexport_pack_label(string scalar label, real scalar j)
{
    real colvector values, lengths
    string colvector text, entries
    string scalar packed, stub
    real scalar width

    st_local("__cx_ok", "0")
    if (!st_vlexists(label)) return
    st_vlload(label, values, text)
    if (rows(values) == 0) return
    lengths = strlen(text)
    width = max((1 \ lengths))
    if (width > 2045) return
    entries = strtrim(strofreal(values, "%21.0f")) :+ " " :+ strofreal(lengths) :+ " " :+ text
    packed = strofreal(rows(values)) + " " + invtokens(entries', " ") + " |"
    if (strlen(packed) != strlen(strofreal(rows(values))) + sum(strlen(entries)) + rows(values) + 2) return
    if (strlen(packed) >= c("macrolen") | strlen(packed) > 1000000) return
    stub = "__cexport_map_" + strofreal(j)
    st_local(stub, packed)
    st_local(stub + "_len", strofreal(strlen(packed), "%21.0f"))
    st_local("__cexport_mapwidth_" + strofreal(j), strofreal(width))
    st_local("__cx_ok", "1")
}
end

* The programs below were separate ado files until ctools 1.3.2. Stata's
* net install takes at most 100 files per package, so helpers that only
* this file calls are kept here as its subprograms.

* Binary statistical writers; serialization is provided by C.
program define _cio_export, rclass
    version 14.1
    gettoken kind 0 : 0
    local syntax_opts "REPLACE"
    if "`kind'"=="spss" local syntax_opts "`syntax_opts' NOVALLabels"
    * export sasxport5 takes val() for vallabfile()
    if "`kind'"=="sasxport5" local syntax_opts "`syntax_opts' REName VALlabfile(string)"
    if "`kind'"=="sasxport8" local syntax_opts "`syntax_opts' VALLabfile"
    if "`kind'"=="dbase" local syntax_opts "`syntax_opts' VERsion(string) DATAFmt ORIGDBFDATE"
    capture syntax [varlist] [if] [in] using/ [, `syntax_opts']
    if _rc {
        local 0 `"using `macval(0)'"'
        syntax [if] [in] using/ [, `syntax_opts']
    }
    * Stata answers an if one observation at a time. Without one the plugin
    * takes every observation in range; an if is marked and passed last, as a
    * sample marker the plugin reads in parallel. The plugin reports the count.
    local touse ""
    local __cio_select "skipif"
    if `"`if'"' != "" {
        marksample touse, novarlist
        local __cio_select "touse"
    }
    if "`varlist'" == "" unab varlist : _all
    * dBase limits, checked as export dbase does: the variable count of the
    * selection and of the dataset, names, string widths, then the
    * observation width of the selection (c(width) after keeping it).
    if "`kind'"=="dbase" {
        if `: word count `varlist'' > 255 | c(k) > 255 {
            di as error "too many variables specified"
            di as error "    Number of variables must be less than 255."
            exit 103
        }
        local width 0
        foreach v of local varlist {
            if strlen("`v'") > 10 {
                di as error "`v': variable name must be less than 10 characters"
                exit 693
            }
            local type : type `v'
            if "`type'" == "strL" {
                di as error "variable `v' cannot be a strL"
                exit 693
            }
            if substr("`type'",1,3) == "str" & real(substr("`type'",4,.)) > 254 {
                di as error "variable `v' cannot exceed 254 characters in length"
                exit 108
            }
            local width = `width' + cond(substr("`type'",1,3)=="str", real(substr("`type'",4,.)), cond("`type'"=="byte",1,cond("`type'"=="int",2,cond("`type'"=="double",8,4))))
        }
        if `width' > 4000 {
            di as error "invalid observation width"
            di as error "    Observation width must be less than 4000."
            exit 902
        }
    }
    local ext ".sas7bdat"
    if "`kind'" == "dbase" local ext ".dbf"
    if "`kind'" == "spss" local ext ".sav"
    if "`kind'" == "sasxport5" local ext ".xpt"
    if "`kind'" == "sasxport8" local ext ".v8xpt"
    * ~ is the home directory, as for Stata's own file commands
    _ctools_home using `"`macval(using)'"'
    mata: st_local("base",pathbasename(st_local("using")))
    if !strpos(`"`macval(base)'"',".") {
        local using `"`macval(using)'`ext'"'
    }
    if "`replace'" == "" {
        confirm new file `"`macval(using)'"'
    }
    * Native's only explicit legacy version is version(3); IV is the default.
    if "`version'" != "" & "`version'" != "3" exit 198
    if "`version'" == "3" local version "III"
    local __cio_rename = ("`rename'" != "")
    local __cio_version "`version'"
    local __cio_origdbfdate = ("`origdbfdate'" != "")
    foreach field in year month day {
        local __cio_dbf_`field' : char _dta[dbf_`field']
    }
    local __cio_datafmt = ("`datafmt'" != "")
    local __cio_xpfoutput ""
    if "`kind'"=="sasxport5" {
        if "`vallabfile'"=="" local vallabfile "xpf"
        foreach full in xpf sascode both none {
            if "`vallabfile'" != "" & substr("`full'",1,strlen("`vallabfile'")) == "`vallabfile'" local vallabfile "`full'"
        }
        if !inlist("`vallabfile'","xpf","sascode","both","none") exit 198
        if inlist("`vallabfile'","xpf","both") {
            * formats.xpf beside the file, using only Stata 14.1 Mata functions.
            mata: st_local("__cio_xpfoutput",substr(st_local("using"),1,strlen(st_local("using"))-strlen(pathbasename(st_local("using"))))+"formats.xpf")
        }
    }
    local __cio_sascode ""
    if ("`kind'"=="sasxport8" & "`vallabfile'"!="") | ("`kind'"=="sasxport5" & inlist("`vallabfile'","sascode","both")) {
        mata: st_local("__cio_sascode",pathrmsuffix(st_local("using"))+".sas")
    }
    mata: st_local("__cio_cwd",pwd())
    mata: st_local("__cio_transportpath",pathisabs(st_local("using"))?st_local("using"):pathjoin(pwd(),st_local("using")))
    local __cio_filename : copy local using
    local __cio_replace = ("`replace'" != "")
    local __cio_datalabel : data label
    mata: st_local("__cio_table",strupper(pathrmsuffix(pathbasename(st_local("using")))))
    * Native names a version 5 member after the file name up to its first
    * period, unchanged otherwise (the plugin keeps eight bytes).
    if "`kind'"=="sasxport5" {
        mata: st_local("__cio_table",strupper(pathbasename(st_local("using"))))
        mata: st_local("__cio_table",substr(st_local("__cio_table"),1,strpos(st_local("__cio_table")+".",".")-1))
    }
    local j=0
    foreach v of local varlist {
        local ++j
        local __cio_name_`j' "`v'"
        local __cio_type_`j' : type `v'
        local __cio_format_`j' : format `v'
        local __cio_label_`j' : variable label `v'
        local __cio_labelset_`j' ""
        local __cio_lcount_`j' 0
        if "`novallabels'" == "" {
            local __cio_labelset_`j' : value label `v'
            if "`__cio_labelset_`j''" != "" {
                mata: _cio_label_metadata(`j',st_local("__cio_labelset_`j'"))
            }
        }
    }
    local __cio_nvar `j'
    local referenced ""
    forvalues i=1/`j' {
        local referenced `referenced' `__cio_labelset_`i''
    }
    quietly label dir
    local allsets `r(names)'
    * Version 5 transport writes, and checks the names of, only the value
    * labels that variables use, as export sasxport5 does.
    if "`kind'"=="sasxport8" {
        foreach set of local allsets {
            if !`: list set in referenced' {
                local ++j
                local __cio_labelset_`j' "`set'"
                mata: _cio_label_metadata(`j',st_local("set"))
            }
        }
    }
    local __cio_nmeta `j' 
    local exportvars `varlist'
    * SPSS and SAS version 8 transport timestamps exclude leap seconds, and
    * both store monthly, quarterly, weekly and half-yearly dates as the date
    * of the period's first day. Keep source metadata and convert only the
    * values passed to the writer; version 5 transport retains Stata units.
    if inlist("`kind'","spss","sasxport8") {
        local exportvars ""
        foreach v of local varlist {
            local fmt : format `v'
            local unit ""
            if regexm("`fmt'","^%-?t([Cmqwh])") local unit = regexs(1)
            if "`unit'" != "" {
                local convert = cond("`unit'"=="C","cofC",cond("`unit'"=="m","dofm",cond("`unit'"=="q","dofq",cond("`unit'"=="w","dofw","dofh"))))
                tempvar ordinary_time
                quietly gen double `ordinary_time' = `convert'(`v') if !missing(`v')
                quietly replace `ordinary_time' = `v' if missing(`v')
                local exportvars `exportvars' `ordinary_time'
            }
            else local exportvars `exportvars' `v'
        }
    }
    if "`kind'"=="dbase" & "`datafmt'"!="" {
        local j=0
        local formatted ""
        foreach v of local varlist {
            local ++j
            local fmt : format `v'
            if !strpos("`: type `v''","str") & substr("`fmt'",1,2)!="%t" & "`fmt'"!="%d" {
                if !regexm("`fmt'","^%-?[0-9]+[.][0-9]+[fge]$") exit 120
                tempvar dbffmt`j'
                quietly gen str2045 `dbffmt`j'' = strtrim(string(`v',"`fmt'"))
                local formatted `formatted' `dbffmt`j''
                local __cio_fmtvar_`j' = `__cio_nvar' + `: word count `formatted''
            }
        }
        local exportvars `exportvars' `formatted'
    }
    _ctools_load
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    capture noisily plugin call ctools_plugin `exportvars' `touse' `in', "cio export `kind' `__cio_select'"
    * A write-protected destination (plugin r(608)) gives native's code and text
    if _rc == 608 {
        if inlist("`kind'","spss","sasxport8") {
            di as error `"file `macval(using)' already exists"'
        }
        else if "`kind'"=="dbase" {
            di as error `"file `macval(using)' could not be opened"'
        }
        else {
            di as error `"file `macval(using)' is read-only; cannot be modified or erased"'
        }
        exit cond(inlist("`kind'","spss","sasxport8","dbase"),602,608)
    }
    if _rc exit _rc
    return local filename : copy local using
    return scalar N = `__cio_N'
    return scalar k = `__cio_nvar'
    di as text `"file `macval(using)' saved"'
end
mata:
void _cio_label_metadata(real scalar j,string scalar set) {
    real colvector values
    string colvector text
    real scalar k
    st_vlload(set,values,text)
    st_local(sprintf("__cio_lcount_%g",j),strofreal(rows(values)))
    for (k=1;k<=rows(values);k++) {
        st_local(sprintf("__cio_lv_%g_%g",j,k),strofreal(values[k],"%24.17e"))
        st_local(sprintf("__cio_lt_%g_%g",j,k),text[k])
    }
}
end

program define _cio_shp_export, rclass
    version 14.1
    syntax [anything] [using/] [, REPLACE SHX ID(varname numeric)]
    if `"`using'"'=="" gettoken using extra : anything
    if `"`using'"'=="" exit 198
    mata: st_local("base",pathbasename(st_local("using")))
    if !strpos(`"`base'"',".") local using `"`using'.shp"'
    if "`id'"=="" local id _ID
    confirm numeric variable `id' _X _Y
    confirm string variable rec_header
    * The characteristic may come from someone else's file: read it in Mata
    * and accept only a number, never re-expanding its text.
    mata: st_local("shape",strtrim(st_global("_dta[shp_type]")))
    mata: st_local("shape",st_local("shape")=="" ? "" : strofreal(strtoreal(st_local("shape")),"%12.0g"))
    if "`shape'"=="" | "`shape'"=="." exit 198
    local __cio_shp_type "`shape'"
    foreach field in xmin ymin xmax ymax zmin zmax mmin mmax {
        local __cio_shp_`field' : char _dta[shp_`field']
    }
    local __cio_filename `"`using'"'
    local __cio_replace = ("`replace'"!="")
    local __cio_shx = ("`shx'"!="")
    _ctools_load
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc!=0 & _rc!=110 exit 601
    local variables `id' _X _Y rec_header
    if inlist(`shape',11,13,15,18,31) local variables `variables' _Z
    if inlist(`shape',11,13,15,18,21,23,25,28,31) local variables `variables' _M
    capture noisily plugin call ctools_plugin `variables', "cio export shp"
    * A write-protected destination (plugin r(608)) gives native's code and text
    if _rc == 608 {
        di as error "could not copy shp file"
        exit 602
    }
    if _rc exit _rc
    return local filename `"`using'"'
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

* ctools: a file name with a leading ~ read as the home directory
* _ctools_home lname name sets the caller's local lname to name, with a
* leading ~ (alone, or before / or \) replaced by the home directory: HOME,
* or USERPROFILE on Windows, as Stata's own file commands read it. The
* plugin opens names as they are, so ~/data.csv would name a folder ~ in the
* working directory. Other names, URLs included, are left as they are.
program define _ctools_home
    version 14.1
    args lname name
    local home : environment HOME
    if c(os) == "Windows" local home : environment USERPROFILE
    mata: st_local("name", _ctools_home_name(st_local("name"), st_local("home")))
    c_local `lname' `"`macval(name)'"'
end

version 14.1
mata:
string scalar _ctools_home_name(string scalar name, string scalar home)
{
    if (home == "" | substr(name, 1, 1) != "~") return(name)
    if (strlen(name) > 1 & !anyof(("/", "\"), substr(name, 2, 1))) return(name)
    return(home + substr(name, 2, .))
}
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
