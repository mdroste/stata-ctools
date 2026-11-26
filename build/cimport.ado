*! version 1.3.2 20261124 github.com/mdroste/stata-ctools
*! original: Stata's import delimited (import_delimited.ado 1.3.5, April 18, 2022)
*! original: Stata's import excel (import_excel.ado 1.0.15, September 20, 2019)
*! original: see help cimport for details

program define cimport, rclass
    version 14.1

    * Parse the subcommand
    gettoken subcmd 0 : 0, parse(" ,")
    if inlist("`subcmd'", "delim", "delimi", "delimit", "delimite") local subcmd "delimited"
    if inlist("`subcmd'", "exc", "exce") local subcmd "excel"

    if "`subcmd'" == "sasxport" & _caller()<16 local subcmd "sasxport5"
    if inlist("`subcmd'","hav","have") local subcmd "haver"
    if inlist("`subcmd'","haverd","haverdi","haverdir","haverdire","haverdirec") local subcmd "haverdirect"

    * FRED (a web API) and Haver (Haver Analytics' Windows DLX library) are
    * online services, not file formats: Stata's own commands handle them,
    * with every argument passed through unchanged. Options snapshot and fast
    * do not apply; the native command manages the data in memory itself.
    if inlist("`subcmd'","fred","haver","haverdirect") {
        version `=_caller()': import `subcmd' `macval(0)'
        return add
        exit
    }
    * macval(): text such as sheet("Q1\$Rev") reaches the handler as typed,
    * not with its $ or backquotes expanded a second time.
    if inlist("`subcmd'", "sas", "spss", "sasxport5", "sasxport8", "dbase", "shp") {
        _cio_import `subcmd' `macval(0)'
        return add
        exit
    }

    if "`subcmd'" == "parquet" {
        _cimport_parquet `macval(0)'
        return add
        exit
    }

    if "`subcmd'" == "excel" {
        * Dispatch to Excel import handler
        cimport_excel `macval(0)'
        return add
        exit
    }
    else if "`subcmd'" != "delimited" {
        di as error "cimport: unknown subcommand `subcmd'"
        di as error "Supported file formats: delimited, excel, parquet, sas, spss, sasxport5, sasxport8, dbase, shp"
        exit 198
    }

    * Now parse the rest with import delimited-style syntax
    * Support both "using filename" and just "filename" (like import delimited does)
    syntax [anything] [using/] [, Delimiters(string asis) VARNames(string) CLEAR ///
        CASE(string) ENCoding(string) BINDQuotes(string) ///
        STRIPQuotes(string) ROWRange(string) COLRange(string) Verbose THReads(integer 0) ///
        ASFloat ASDOUBle NUMERICcols(string) STRINGcols(string) ///
        DECIMALSeparator(string) GROUPSeparator(string) ///
        MAXQUOTEDrows(string) EMPTYlines(string) ///
        PARSELocale(string) CHARSET(string) FAVORSTRFixed COLLAPSEDelimiters ///
        SNAPshot FAST TEMPDTA]
    if "`snapshot'" != "" & "`fast'" != "" {
        di as error "options snapshot and fast may not be combined"
        exit 198
    }

    if inlist("`stripquotes'","y","ye") local stripquotes "yes"
    if "`stripquotes'"=="n" local stripquotes "no"
    if strlen("`stripquotes'")>=3 & substr("default",1,strlen("`stripquotes'"))=="`stripquotes'" local stripquotes "default"
    if !inlist("`stripquotes'", "", "default", "yes", "no") {
        di as error "stripquotes() must be default, yes, or no"
        exit 198
    }
    local __cimport_stripquotes "`stripquotes'"
    local __cimport_favorstrfixed = ("`favorstrfixed'"!="")
    if "`maxquotedrows'" == "" local maxquotedrows 20
    if "`maxquotedrows'" == "unlimited" local maxquotedrows 0
    capture confirm integer number `maxquotedrows'
    if _rc | `"`maxquotedrows'"' == "" {
        di as error "option {bf:maxquotedrows()} incorrectly specified."
        exit 198
    }
    * Native accepts any integer: a negative one, like 0, sets no limit
    if `maxquotedrows'<0 local maxquotedrows 0
    foreach opt in numericcols stringcols {
        if "``opt''" != "" & "``opt''" != "_all" {
            capture numlist "``opt''", integer range(>0) sort
            if _rc {
                di as error "{bf:``opt''}: invalid numlist in {bf:`opt'()} option"
                exit 198
            }
            local `opt' "`r(numlist)'"
        }
    }
    if "`numericcols'" == "_all" & "`stringcols'" == "_all" {
        di as error "numericcols(_all) and stringcols(_all) cannot be combined"
        exit 198
    }

    * Handle filename - can come from using/ or as first positional argument
    local explicit_names ""
    if `"`using'"' != "" & `"`anything'"' != "" {
        local explicit_names `"`anything'"'
    }
    if `"`using'"' == "" & `"`anything'"' != "" {
        * Filename provided without "using" keyword
        gettoken using extra : anything
        if trim(`"`extra'"') != "" {
            di as error "cimport: specify one filename"
            exit 198
        }
    }
    if `"`using'"' == "" {
        di as error "cimport delimited: filename required"
        di as error "Syntax: cimport delimited [using] filename [, options]"
        exit 198
    }

    * As import delimited does, a URL is copied to a temporary file first, as
    * it is: no .csv is added to it, and r(filename) reports the URL.
    local __cimport_source `"`macval(using)'"'
    mata: st_local("__isurl", strofreal(pathisurl(st_local("using"))))
    if `__isurl' {
        tempfile __urlcopy
        capture copy `"`macval(using)'"' `"`__urlcopy'"'
        if _rc {
            di as error "could not open url"
            exit 603
        }
        local using `"`__urlcopy'"'
    }
    else {
        * ~ is the home directory, as for Stata's own file commands
        _ctools_home using `"`macval(using)'"'
        * Native delimited import supplies .csv when no extension is given.
        mata: st_local("__basename", pathbasename(st_local("using")))
        if !strpos(`"`macval(__basename)'"', ".") {
            local using `"`macval(using)'.csv"'
        }
        local __cimport_source `"`macval(using)'"'
    }
    * Validate file exists
    confirm file `"`macval(using)'"'

    * Like native, explicit names must be valid; they rename variables
    * 1..k after the import, so header detection stays automatic.
    foreach name of local explicit_names {
        mata: st_local("valid", strofreal(st_isvarname(st_local("name"))))
        if !`valid' {
            di as error "{bf:`name'}: invalid variable name"
            exit 198
        }
    }

    * Scan and validate before discarding the current dataset.
    if "`clear'" == "" & (_N > 0 | c(k) > 0) {
        di as error "data in memory would be lost"
        di as error "use the clear option to discard current data"
        exit 4
    }

    local __cimport_autoheader = ("`varnames'" == "")

    * Set default delimiter (auto-detect if not specified, matching Stata behavior)
    local plugin_delim "auto"
    if `"`delimiters'"' != "" {
        local plugin_delim `"`delimiters'"'
    }

    local __cimport_delimiters ""
    local __cimport_collapse 0
    local __cimport_asstring 0
    if `"`delimiters'"' != "" {
        _cimport_delimiters `macval(delimiters)'
        local delimiters `"`r(delimiters)'"'
        local __cimport_delimiters `"`r(delimiters)'"'
        local __cimport_collapse = r(collapse)
        local __cimport_asstring = r(asstring)
        if strlen(`"`delimiters'"') == 1 & !`__cimport_collapse' {
            local plugin_delim `"`delimiters'"'
            if `"`delimiters'"' == char(9) local plugin_delim "tab"
            if `"`delimiters'"' == " " local plugin_delim "space"
        }
        else local plugin_delim "auto"
    }
    if "`charset'" != "" {
        if "`encoding'" == "" local encoding "`charset'"
    }

    if strlen("`varnames'")>=3 & substr("nonames",1,strlen("`varnames'"))=="`varnames'" local varnames "nonames"
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
            if _rc != 0 | real("`varnames'") < 1 {
                di as error "cimport: varnames() must be a positive integer or nonames"
                exit 198
            }
            local headerrow = `varnames'
        }
    }

    if "`case'"!="" {
        if substr("lower",1,strlen("`case'"))=="`case'" local case "lower"
        if substr("upper",1,strlen("`case'"))=="`case'" local case "upper"
        if strlen("`case'")>=3 & substr("preserve",1,strlen("`case'"))=="`case'" local case "preserve"
    }
    * Parse case option (default is lower)
    if "`case'" != "" {
        if !inlist("`case'", "preserve", "lower", "upper") {
            di as error "cimport: case() must be preserve, lower, or upper"
            exit 198
        }
    }
    else {
        local case "lower"
    }

    * The C converter validates the supplied charset and never silently guesses.
    local encoding=strtrim(`"`macval(encoding)'"')
    local __cimport_encoding_name `"`encoding'"'
    local encoding_opt ""

    if "`bindquotes'"!="" {
        foreach full in loose strict nobind {
            if substr("`full'",1,strlen("`bindquotes'"))=="`bindquotes'" local bindquotes "`full'"
        }
    }
    * Parse bindquotes - default is loose (matches Stata's import delimited default)
    if "`bindquotes'" != "" {
        if !inlist("`bindquotes'", "strict", "loose", "nobind") {
            di as error "cimport: bindquotes() must be strict, loose, or nobind"
            exit 198
        }
    }
    else {
        local bindquotes "loose"
    }

    * Parse rowrange option
    * Native rowrange() counts 1-based file lines (the header is a line too).
    * The lines go to C unchanged; it maps them to rows once it has decided
    * whether there is a header.
    local startrow = 0
    local endrow = 0
    if "`rowrange'" != "" {
        _cimport_range `"`rowrange'"' row
        local startrow = r(first)
        local endrow = r(last)
    }

    local __cimport_startrow `startrow'
    local __cimport_endrow `endrow'
    * Parse colrange option
    local startcol = 0
    local endcol = 0
    if "`colrange'" != "" {
        _cimport_range `"`colrange'"' column
        local startcol = r(first)
        local endcol = r(last)
    }
    * C applies colrange() before header detection and naming, as native does.
    local __cimport_startcol `startcol'
    local __cimport_endcol `endcol'

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

    * emptylines(): native accepts exactly skip or include. Its reader follows
    * the caller's version: from version 17, include keeps each empty line
    * as an observation (the plugin's emptylines=fill); before 17 the older
    * reader reads include as skip.
    if !inlist(`"`emptylines'"', "", "skip", "include") {
        di as error `"{bf:`emptylines'} is an invalid {bf:emptylines()} argument; use {bf:skip}, or {bf:include}"'
        exit 198
    }
    if "`emptylines'" == "include" & _caller() >= 17 local emptylines "fill"
    else local emptylines "skip"

    * Numeric locale profiles and Unicode digits are interpreted entirely in C.
    local __cimport_parselocale "`parselocale'"
    local __cimport_decimalseparator `"`decimalseparator'"'
    local __cimport_groupseparator `"`groupseparator'"'

    * Validate maxquotedrows
    if `maxquotedrows' < 0 {
        di as error "cimport: maxquotedrows() must be non-negative"
        exit 198
    }

    * numericcols()/stringcols() reach C as space-separated column numbers
    * (or _all) with their byte lengths, so a long list is never cut short
    local __cimport_numcols "`numericcols'"
    local __cimport_strcols "`stringcols'"
    mata: st_local("__cimport_numcols_len", strofreal(strlen(st_local("numericcols"))))
    mata: st_local("__cimport_strcols_len", strofreal(strlen(st_local("stringcols"))))

    * Load the platform-appropriate ctools plugin if not already loaded
    _ctools_load
    * Stata scopes plugin registrations to the calling ado program.
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    capture confirm number 0

    * Build plugin arguments
    local opt_noheader = cond(`noheader' == 1, "noheader", "")
    local opt_headerrow = cond(`headerrow' > 1, "headerrow=`headerrow'", "")
    local opt_verbose = cond("`verbose'" != "", "verbose", "")
    local opt_case = "case=`case'"
    local opt_bindquotes = "bindquotes=`bindquotes'"

    * New options
    local opt_asfloat = cond("`asfloat'" != "", "asfloat", "")
    local opt_asdouble = cond("`asdouble'" != "" | ("`asfloat'" == "" & "`c(type)'" == "double"), "asdouble", "")
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

    * Build threads option string
    local threads_code ""
    if `threads' > 0 {
        local threads_code "threads(`threads')"
    }

    * Header and transfer globals of an earlier, interrupted import
    capture macro drop _cimport_label_* _cimport_blob*
    * Timing scalars come only from a load: none may be left by an earlier run
    capture scalar drop _cimport_time_mmap _cimport_time_encoding _cimport_time_parse _cimport_time_infer
    capture scalar drop _cimport_time_cache _cimport_time_store _cimport_time_total
    capture scalar drop _cimport_threads_max _cimport_openmp_enabled

    * Record start time (the plugin is loaded, so its clock is called directly)
    plugin call ctools_plugin, "clock"
    local __ct_start99 `ctools_clock_value'


    * =========================================================================
    * PHASE 1: Scan the CSV to get metadata
    * =========================================================================

    local __cimport_filename `"`macval(using)'"'
    capture noisily plugin call ctools_plugin, ///
        "cimport `threads_code' scan @filename `plugin_delim' `opt_noheader' `opt_headerrow' `opt_verbose' `opt_bindquotes' `opt_asfloat' `opt_asdouble' `opt_decimalsep' `opt_groupsep' `opt_emptylines' `opt_maxquotedrows' `encoding_opt'"

    local scan_rc = _rc
    if `scan_rc' {
        di as error "Error scanning CSV file (rc=`scan_rc')"
        capture macro drop _cimport_label_*
        exit `scan_rc'
    }

    * Retrieve metadata from global macros (SF_macro_save creates globals)
    local nobs = ${_cimport_nobs}
    local nvar = ${_cimport_nvar}
    local varnames ${_cimport_varnames}
    local vartypes ${_cimport_vartypes}
    local numtypes ${_cimport_numtypes}
    local strlens ${_cimport_strlens}
    local longtypes ${_cimport_longtypes}
    local strlbytes ${_cimport_strlbytes}
    local hasheader = 0${_cimport_hasheader}
    mata: st_local("returned_encoding",st_global("_cimport_encoding"))
    mata: st_local("returned_delimiters",subinstr(st_global("_cimport_delimiters"),char(9),char(92)+"t"))

    * Clean up global macros
    macro drop _cimport_nobs _cimport_nvar _cimport_varnames ///
               _cimport_vartypes _cimport_numtypes _cimport_strlens _cimport_longtypes _cimport_strlbytes ///
               _cimport_encoding _cimport_delimiters _cimport_hasheader

    * Like import delimited, name the charset only when it was detected.
    if `"`__cimport_encoding_name'"' == "" di as text `"(encoding automatically selected: `macval(returned_encoding)')"'

    * From here every exit before the load releases the scanned file.
    * More columns than maxvar allows fail as native does (after colrange()).
    if `nvar' > c(maxvar) {
        di as error "There are more columns than allowed by this Stata. See help maxvar for more information."
        capture plugin call ctools_plugin, "cimport clear"
        capture macro drop _cimport_label_*
        exit 900
    }
    * More explicit names than columns fails before loading (rechecked after colrange()).
    if `: word count `explicit_names'' > `nvar' {
        di as error "too many variables specified"
        capture plugin call ctools_plugin, "cimport clear"
        capture macro drop _cimport_label_*
        exit 103
    }
    if `nobs' > 2147483647 {
        di as error "ctools does not support datasets exceeding 2^31 (2.147 billion) observations"
        di as error "This is a limitation of Stata's plugin API"
        capture plugin call ctools_plugin, "cimport clear"
        capture macro drop _cimport_label_*
        exit 920
    }
    * A failed or interrupted import leaves the data in memory as they were.
    * The import replaces the whole dataset, so the old data are copied first:
    * from Stata 16 into a frame of cimport's own, in memory (the copy fails
    * before anything changes if memory is short), and restored from it with
    * frame copy. The import is not built in a new frame that then takes the
    * current frame's name: when a frame holding native import excel data is
    * renamed or dropped and another frame takes its name, a later native
    * import spss returns wrong strings (StataNow 19.5).
    * Option snapshot preserves the data instead (on disk above
    * c(max_preservemem)), as do imports before Stata 16 and imports into an
    * empty dataset, whose value labels and characteristics stay without
    * clear. Option fast keeps no copy and undoes nothing. preserve and the
    * frame belong to this program, so these steps are inline.
    local __cimport_frame ""
    local __cimport_snap 0
    * The strL transfer makes a frame of its own current while it copies a
    * block; a Break there must not restore into that frame.
    local __cimport_home ""
    if c(stata_version) >= 16 local __cimport_home = c(frame)
    if "`fast'" == "" {
        if "`snapshot'" == "" & c(stata_version) >= 16 & (_N > 0 | c(k) > 0) {
            tempname __cimport_frame
            capture noisily frame copy `c(frame)' `__cimport_frame'
        }
        else {
            capture noisily preserve
            local __cimport_snap 1
        }
        if _rc {
            local rc = _rc
            capture plugin call ctools_plugin, "cimport clear"
            capture macro drop _cimport_label_*
            exit `rc'
        }
    }
    if "`clear'" != "" clear

    * Handle empty file or file with header only - match Stata's behavior (rc=0, N=0, k=0)
    if `nvar' == 0 {
        capture plugin call ctools_plugin, "cimport clear"
        plugin call ctools_plugin, "clock"
        local elapsed = `ctools_clock_value' - `__ct_start99'
        di as text "(0 vars, 0 obs)"
        return scalar N = 0
        return scalar k = 0
        return scalar time = `elapsed'
        return local filename `"`macval(__cimport_source)'"'
        return local encoding `"`macval(returned_encoding)'"'
        return local delimiters `"`macval(returned_delimiters)'"'
        if "`__cimport_frame'" != "" frame drop `__cimport_frame'
        else if `__cimport_snap' restore, not
        exit 0
    }
    * Phases 2 and 3 run as one block: a failure or Break in them restores
    * the data in memory
    capture noisily {


    * =========================================================================
    * PHASE 2: Create variables in Stata
    * =========================================================================

    * Names, types, labels and formats are settled in Mata and all variables
    * are created in one call; strL columns are returned in local longcolumns.
    * With tempdta, when the import may replace the dataset outright, with
    * nothing in memory that use would discard (value labels, a data label or
    * dataset characteristics of an empty dataset kept without clear), and
    * strL text dominates the data, the plugin writes a temporary .dta that
    * use loads: the whole dataset (cimport fulldta), or failing that its strL
    * columns (cimport strldta), several times faster than Mata can store
    * strLs. tempdta is an option because the file, about the size of the
    * imported data, may reach the disk. c(filename) and c(filedate) stay as
    * they were.
    local longcolumns ""
    local __cimport_dta 0
    local __cimport_full 0
    local __cimport_loaded 0
    local __cimport_viadta 0
    if `nobs' & "`tempdta'" != "" {
        quietly label dir
        local __cimport_keep `"`r(names)'"'
        local __cimport_label : data label
        local __cimport_chars : char _dta[]
        if `"`__cimport_keep'`__cimport_label'`__cimport_chars'"' == "" local __cimport_viadta 1
    }
    mata: _cimport_addvars(st_local("case"), `__cimport_viadta')
    if `__cimport_full' {
        tempfile __cimport_blobpath
        mata: st_local("__cimport_fullnames_len", strofreal(strlen(st_local("__cimport_fullnames"))))
        mata: st_local("__cimport_fulltypes_len", strofreal(strlen(st_local("__cimport_fulltypes"))))
        * A refusal or a failed write leaves the scan for the next route.
        * Validation hook: fulldta refuses, as it does when it cannot
        * group the rows.
        if "$CTOOLS_cimport_fault" == "fulldta" local __cimport_fullrc 459
        else {
            plugin call ctools_plugin, "cimport `threads_code' fulldta"
            local __cimport_fullrc = 0${_cimport_fullrc}
            capture macro drop _cimport_fullrc
        }
        if `__cimport_fullrc' == 1 exit 1
        if !`__cimport_fullrc' {
            mata: st_local("__cimport_fn", st_global("S_FN")); st_local("__cimport_fndate", st_global("S_FNDATE"))
            quietly use `"`__cimport_blobpath'"', clear
            mata: st_global("S_FN", st_local("__cimport_fn")); st_global("S_FNDATE", st_local("__cimport_fndate"))
            capture erase `"`__cimport_blobpath'"'
            plugin call ctools_plugin, "cimport clear"
            mata: _cimport_addvars(st_local("case"), 3)
            local __cimport_loaded 1
            local __cimport_dta 0
        }
        else {
            capture erase `"`__cimport_blobpath'"'
            * Nothing was created: take the next route.
            if !`__cimport_dta' mata: _cimport_addvars(st_local("case"), 0)
        }
    }
    if `__cimport_dta' {
        tempfile __cimport_blobpath
        local __cimport_longcolumns "`longcolumns'"
        mata: st_local("__cimport_longcolumns_len", strofreal(strlen(st_local("longcolumns"))))
        mata: st_local("__cimport_strlnames_len", strofreal(strlen(st_local("__cimport_strlnames"))))
        plugin call ctools_plugin, "cimport `threads_code' strldta"
        mata: st_local("__cimport_fn", st_global("S_FN")); st_local("__cimport_fndate", st_global("S_FNDATE"))
        quietly use `"`__cimport_blobpath'"', clear
        mata: st_global("S_FN", st_local("__cimport_fn")); st_global("S_FNDATE", st_local("__cimport_fndate"))
        capture erase `"`__cimport_blobpath'"'
        mata: _cimport_addvars(st_local("case"), 2)
    }
    else if `nobs' & !`__cimport_loaded' quietly set obs `nobs'

    * =========================================================================
    * PHASE 3: Load data into variables
    * =========================================================================

    * strL columns arrive a block of rows at a time and are stored in Mata.
    * Values found in the source, verbatim or with doubled quotes or CR LF to
    * undo, travel as positions, which Mata cuts from the source itself. From
    * Stata 16 the positions come through a frame, with nothing written to a
    * file (see _cimport_read_blobmem); otherwise, or when the source cannot
    * be cut in place, through a short-lived file of blocks.
    if `nobs' & "`longcolumns'" != "" & !`__cimport_dta' & !`__cimport_loaded' {
        local __cimport_longcolumns "`longcolumns'"
        mata: st_local("__cimport_longcolumns_len", strofreal(strlen(st_local("longcolumns"))))
        mata: st_local("__cimport_blobdirect", strofreal(_cimport_canopen(st_local("__cimport_filename"))))
        local __cimport_memrc 459
        * (Validation hook: blobmem takes the file transfer.)
        if c(stata_version) >= 16 & `__cimport_blobdirect' & "$CTOOLS_cimport_fault" != "blobmem" {
            plugin call ctools_plugin, "cimport `threads_code' blobmem"
            local __cimport_memrc = 0${_cimport_blobmemrc}
            if !`__cimport_memrc' mata: _cimport_read_blobmem(st_local("__cimport_filename"), strtoreal(tokens(st_local("longcolumns"))), `nobs')
            capture macro drop _cimport_blobmemrc _cimport_blobblocks _cimport_blobmax _cimport_blobsource ///
                _cimport_blobinfo _cimport_blobcounts _cimport_blobtext _cimport_blobnext
        }
        if `__cimport_memrc' {
            tempfile __cimport_blobpath
            plugin call ctools_plugin, "cimport `threads_code' blobs"
            mata: _cimport_read_blobs(st_local("__cimport_blobpath"), strtoreal(tokens(st_local("longcolumns"))), `nobs')
            capture erase `"`__cimport_blobpath'"'
        }
    }

    if !`nobs' plugin call ctools_plugin, "cimport clear"
    if `nobs' & !`__cimport_loaded' capture noisily plugin call ctools_plugin *, ///
        "cimport `threads_code' load @filename `plugin_delim' `opt_noheader' `opt_headerrow' `opt_verbose' `opt_bindquotes' `opt_asfloat' `opt_asdouble' `opt_decimalsep' `opt_groupsep' `opt_emptylines' `opt_maxquotedrows' `encoding_opt'"

    local load_rc = cond(`__cimport_loaded', 0, _rc)
    if `load_rc' {
        di as error "Error loading CSV data (rc=`load_rc'); imported data are incomplete"
        exit `load_rc'
    }

    * Native shows a numeric column in the format of the type it infers from
    * a sample of lines (see _cimport_sample_formats), whatever the final
    * storage.
    if `nobs' mata: _cimport_sample_formats(`hasheader', "`opt_asdouble'" != "")

    * Like native, explicit names rename imported variables 1..k as typed
    * (no case conversion; header labels are kept). Errors restore the data.
    if `"`explicit_names'"' != "" {
        if `: word count `explicit_names'' > c(k) {
            di as error "too many variables specified"
            exit 103
        }
        local i = 0
        foreach name of local explicit_names {
            local ++i
            mata: st_local("old", st_varname(`i'))
            capture rename `old' `name'
            if _rc {
                di as error "could not name variables specified"
                exit 198
            }
        }
    }
    * Validation hook: fail as Break would once the import is complete
    if "$CTOOLS_cimport_fault" == "load" error 1
    }
    local rc = _rc
    if `rc' {
        * Release the scanned file first; the restore may need the memory
        capture plugin call ctools_plugin, "cimport clear"
        capture macro drop _cimport_label_* _cimport_blob*
        nobreak {
            if "`__cimport_home'" != "" capture frame change `__cimport_home'
            if "`__cimport_frame'" != "" frame copy `__cimport_frame' `c(frame)', replace
            else if `__cimport_snap' restore
        }
        exit `rc'
    }

    plugin call ctools_plugin, "clock"
    local elapsed = `ctools_clock_value' - `__ct_start99'

    * Display summary (matching import delimited output format)
    di as text "(" strtrim(string(c(k),"%20.0fc")) " var" cond(c(k)==1,"","s") ", " strtrim(string(_N,"%20.0fc")) " obs)"

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
        * Without rows to load (a header-only file) no load ran to time
        capture confirm scalar _cimport_time_total
        if !_rc {
            di as text "  C plugin internals:"
            di as text "    Memory map file:        " as result %8.4f _cimport_time_mmap " sec"
            di as text "    Encoding conversion:    " as result %8.4f _cimport_time_encoding " sec"
            di as text "    Parse CSV:              " as result %8.4f _cimport_time_parse " sec"
            di as text "    Type inference:         " as result %8.4f _cimport_time_infer " sec"
            di as text "    Cache conversion:       " as result %8.4f _cimport_time_cache " sec"
            di as text "    Store to Stata:         " as result %8.4f _cimport_time_store " sec"
            di as text "  {hline 53}"
            di as text "    C plugin total:         " as result %8.4f _cimport_time_total " sec"
            di as text "  {hline 53}"
        }
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
    }
    * Clean up timing scalars (written by every load, shown only with verbose)
    capture scalar drop _cimport_time_mmap _cimport_time_encoding _cimport_time_parse _cimport_time_infer
    capture scalar drop _cimport_time_cache _cimport_time_store _cimport_time_total
    capture scalar drop _cimport_threads_max _cimport_openmp_enabled

    if "`__cimport_frame'" != "" frame drop `__cimport_frame'
    else if `__cimport_snap' restore, not

    * Return results
    return scalar N = _N
    return scalar k = c(k)
    return scalar time = `elapsed'
    return local filename `"`macval(__cimport_source)'"'
    return local encoding `"`macval(returned_encoding)'"'
    return local delimiters `"`macval(returned_delimiters)'"'

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

    * Native permits either allstring or allstring(format), with exclusive syntax.
    local original0 `"`macval(0)'"'
    capture syntax [anything(equalok)] [using/] [, ALLSTRING(string) *]
    local allfmt `"`allstring'"'
    local 0 `"`macval(original0)'"'
    * The same two declarations as native: all and allstri are allstring, and
    * only the full allstring(format) carries a format (all(format) is read
    * as neither, as native reads it).
    syntax [anything(equalok)] [using/] [, SHeet(string) CELLRange(string) FIRSTrow ///
        ALLstring(string) ALLstring CASE(string) CLEAR Verbose DEScribe LOCALE(string) DETAIL ///
        SNAPshot FAST]

    * Native describe is an exclusive workbook-inspection syntax.
    if "`describe'"!="" & (`"`macval(sheet)'`cellrange'`firstrow'`allstring'`macval(allfmt)'`case'`clear'`locale'`detail'"'!="") {
        di as error "option {bf:describe} not allowed"
        exit 198
    }

    if `"`allfmt'"' != "" {
        if "`allstring'" == "allstring" opts_exclusive "allstring allstring()"
        mata: st_local("fmtvalid",strofreal(st_isnumfmt(st_local("allfmt"))))
        if !`fmtvalid' {
            di as error `"`macval(allfmt)': invalid numeric display format in option {bf:allstring()}"'
            exit 198
        }
        local allstring "allstring"
    }

    if substr("`cellrange'",-1,1)==":" {
        di as error "invalid lower right cell range in option {bf:cellrange()}"
        exit 198
    }

    * Handle filename - can come from using/ or as first positional argument
    local extvarlist ""
    if `"`using'"' != "" local extvarlist `"`anything'"'
    if `"`using'"' == "" & `"`anything'"' != "" {
        gettoken using extra : anything
        if trim(`"`extra'"') != "" {
            di as error "cimport: specify one filename"
            exit 198
        }
    }
    if `"`using'"' == "" {
        di as error "cimport excel: filename required"
        di as error "Syntax: cimport excel [using] filename [, options]"
        exit 198
    }
    if "`describe'"!="" & `"`extvarlist'"'!="" {
        di as error "option {bf:describe} not allowed"
        exit 198
    }

    * As native import excel does, a URL is first copied to a temporary file,
    * which keeps the URL's .xls, .xlsx or .xlsm (the reader follows it).
    local source `"`macval(using)'"'
    local urlcopy ""
    mata: st_local("isurl", strofreal(pathisurl(st_local("using"))))
    if `isurl' {
        mata: st_local("urlext", strlower(pathsuffix(st_local("using"))))
        if !inlist("`urlext'", ".xls", ".xlsx", ".xlsm") local urlext ".xlsx"
        tempfile urlbase
        local urlcopy `"`urlbase'`urlext'"'
        capture copy `"`macval(using)'"' `"`urlcopy'"'
        if _rc {
            capture erase `"`urlcopy'"'
            di as error "could not open url"
            exit 603
        }
        local using `"`urlcopy'"'
    }
    else {
        * ~ is the home directory, as for Stata's own file commands
        _ctools_home using `"`macval(using)'"'
        * With no extension, native import tries .xls before .xlsx.
        mata: st_local("base",pathbasename(st_local("using")))
        if !strpos(`"`macval(base)'"',".") {
            capture confirm file `"`macval(using)'.xls"'
            if !_rc {
                local using `"`macval(using)'.xls"'
            }
            else {
                local using `"`macval(using)'.xlsx"'
            }
        }
        local source `"`macval(using)'"'
    }
    capture noisily confirm file `"`macval(using)'"'
    if _rc {
        local rc = _rc
        if `"`urlcopy'"' != "" capture erase `"`urlcopy'"'
        exit `rc'
    }
    if `"`extvarlist'"' != "" & ("`case'" != "" | "`firstrow'" != "") {
        if `"`urlcopy'"' != "" capture erase `"`urlcopy'"'
        * As native: its syntax with a list of variables has no case()
        if "`firstrow'" != "" {
            di as error "option {bf:firstrow} not allowed with extvarlist; specify option {bf:cellrange()} to skip the first row"
        }
        else di as error "option {bf:case()} not allowed"
        exit 198
    }

    * As native, the workbook kind comes from the file's first bytes, not its
    * name: an Office Open XML package (.xlsx, .xlsm, .xltx or any other
    * name; macros are not read) or a binary (BIFF) workbook. Other content
    * keeps the kind its extension names, so the reader reports the damage;
    * with neither, it fails as native does.
    mata: st_local("kind", _cimport_excel_kind(st_local("using")))
    if "`kind'" == "" {
        if lower(substr(`"`macval(using)'"',-4,4)) == ".xls" local kind "xls"
        if inlist(lower(substr(`"`macval(using)'"',-5,5)), ".xlsx", ".xlsm") local kind "xlsx"
    }
    if "`kind'" == "" {
        if `"`urlcopy'"' != "" capture erase `"`urlcopy'"'
        mata: printf("{err}file %s could not be loaded\n", st_local("source"))
        exit 603
    }

    * macval(): a sheet name or file name holding $ or ` reaches the reader as
    * typed. Blocks, not one-line ifs: Stata expands the command of a one-line
    * if a second time, macval() or not.
    local opts ""
    if `"`macval(allfmt)'"' != "" {
        local opts `"allstringfmt(`"`macval(allfmt)'"')"'
    }
    if `"`macval(sheet)'"' != "" {
        local opts `"`macval(opts)' sheet(`"`macval(sheet)'"')"'
    }
    if "`cellrange'" != "" {
        local opts `"`macval(opts)' cellrange(`cellrange')"'
    }
    if "`case'" != "" {
        local opts `"`macval(opts)' case(`case')"'
    }
    capture noisily _cio_import `kind' `extvarlist' using `"`macval(using)'"', ///
        `firstrow' `allstring' `clear' `describe' `snapshot' `fast' `macval(opts)'
    local rc = _rc
    if `"`urlcopy'"' != "" capture erase `"`urlcopy'"'
    if `rc' exit `rc'
    return add
    return local filename `"`macval(source)'"'
    exit

end

program define _cimport_delimiters, rclass
    version 14.1
    * Messages as native import delimited gives them
    capture syntax anything(name=delimiters) [, COLLApse ASSTRing]
    if _rc {
        di as error "invalid {bf:delimiter()} argument"
        exit 198
    }
    if "`collapse'"!="" & "`asstring'"!="" {
        di as error "invalid {bf:delimiter()} argument"
        exit 198
    }
    local quoted = strpos(`"`macval(delimiters)'"',char(34))>0
    gettoken delimiters extra : delimiters
    if strtrim(`"`extra'"')!="" | `"`delimiters'"'=="" {
        di as error "invalid {bf:delimiter()} argument"
        exit 198
    }
    local delimiters = subinstr(`"`macval(delimiters)'"',"\t",char(9),.)
    if !`quoted' {
        if !inlist("`delimiters'","tab","comma","space","whitespace",char(9),",") {
            di as error `"{bf:`delimiters'}: invalid {bf:delimiter()} shortcut"'
            exit 198
        }
        if "`delimiters'"=="tab" local delimiters=char(9)
        if "`delimiters'"=="comma" local delimiters=","
        if "`delimiters'"=="space" local delimiters=" "
        if "`delimiters'"=="whitespace" local delimiters=char(9)+" "
    }
    return local delimiters `"`macval(delimiters)'"'
    return scalar collapse = ("`collapse'"!="")
    return scalar asstring = ("`asstring'"!="")
end

mata:
/* "xlsx" for a zip package (Office Open XML), "xls" for an OLE2 compound
   file (BIFF), "" otherwise or when the file cannot be read. */
string scalar _cimport_excel_kind(string scalar path) {
    real scalar fh
    string scalar head
    fh = _fopen(path, "r")
    if (fh < 0) return("")
    head = fread(fh, 8)
    fclose(fh)
    if (head == J(0, 0, "")) return("")
    if (substr(head, 1, 4) == "PK" + char((3, 4))) return("xlsx")
    if (head == char((208, 207, 17, 224, 161, 177, 26, 225))) return("xls")
    return("")
}
/* cname is the plugin's name for an ASCII header. Other headers keep only
   Unicode letters, decimal digits and underscores, drop leading digits and
   are cut to 32 characters, as native import delimited names them. */
string scalar _cimport_name(string scalar cname, string scalar header) {
    string scalar name
    if (ustrinvalidcnt(cname) == 0 & ustrlen(cname) == strlen(cname)) return(cname)
    name = ustrregexrf(ustrregexra(header, "[^\p{L}\p{Nd}_]", ""), "^\p{Nd}+", "")
    return(usubstr(name, 1, 32))
}
/* Creates the scanned columns as variables in one st_addvar() call. The scan
   leaves the plugin's names and the column types in locals varnames,
   vartypes (1 = string), numtypes (0 double, 1 float, 2 long, 3 int, 4 byte),
   strlens and longtypes (1 = strL), and each header in global
   _cimport_label_#. A name already used, or one that is not a valid variable
   name, becomes v# (V# with case(upper)), then v#_1, v#_2, ...; as in native,
   _all and _skip are valid. A header that differs from its name becomes the
   label. The strL column numbers are left in local longcolumns.

   mode 0 creates every variable. mode 1 (option tempdta) may create nothing
   instead, when there are strL columns and no name is one a varlist reads
   otherwise (_all, _skip), for a temporary .dta that use loads. It leaves 1
   in local __cimport_full, with every name and type in locals
   __cimport_fullnames and __cimport_fulltypes, when at most 32,767
   variables can go in a format-118 dataset (cimport fulldta writes it; mode
   3 then sets formats and labels). It leaves 1 in local __cimport_dta, with
   the strL names in local __cimport_strlnames, when other variables exist
   and a single order command can name them all (cimport strldta writes the
   strL columns; mode 2 then creates the other variables and puts all in file
   order, which rewrites the observations twice). Either route copies the
   other variables' data through a file and use, which costs more than
   storing them through the SPI, so mode 1 takes one only when the strL text
   (local strlbytes, from the scan) is at least as large as the other
   variables' data (local nobs times their width); otherwise Mata storing the
   strLs costs less. */
void _cimport_addvars(string scalar lettercase, real scalar mode)
{
    string rowvector names, vtype, ntype, vlen, vlong, types, formats, labels, longcols
    string scalar name, prefix, label
    real scalar k, i, suffix, width, strl, nostrl, rowwidth, worth
    string scalar tag
    real rowvector index
    transmorphic scalar used

    names = tokens(st_local("varnames"))
    vtype = tokens(st_local("vartypes"))
    ntype = tokens(st_local("numtypes"))
    vlen = tokens(st_local("strlens"))
    vlong = tokens(st_local("longtypes"))
    k = cols(names)
    if (!k) return
    types = formats = labels = J(1, k, "")
    longcols = J(1, 0, "")
    prefix = lettercase == "upper" ? "V" : "v"
    used = asarray_create()
    for (i = 1; i <= k; i++) {
        /* The plugin names ASCII headers. A header with other characters is
           named as native does (see _cimport_name), then the case applies. */
        label = st_global("_cimport_label_" + strofreal(i))
        if (mode != 1) st_global("_cimport_label_" + strofreal(i), "")
        name = _cimport_name(names[i], label)
        if (name == "") name = "v" + strofreal(i)
        if (lettercase == "lower") name = ustrlower(name)
        else if (lettercase == "upper") name = ustrupper(name)
        if (asarray_contains(used, name) | !st_isvarname(name)) {
            name = prefix + strofreal(i)
            for (suffix = 1; asarray_contains(used, name); suffix++) {
                name = prefix + strofreal(i) + "_" + strofreal(suffix)
            }
        }
        asarray(used, name, 1)
        names[i] = name
        /* Header text from the file is stored in Mata only, never pasted
           into a Stata command line. */
        if (label != name) labels[i] = label
        if (vtype[i] == "1") {
            width = strtoreal(vlen[i])
            if (vlong[i] == "1") {
                types[i] = "strL"
                formats[i] = "%9s"
                longcols = longcols, strofreal(i)
            }
            else {
                types[i] = "str" + strofreal(max((1, width)))
                formats[i] = "%" + strofreal(max((9, width))) + "s"
            }
        }
        else if (ntype[i] == "4") {
            types[i] = "byte"
            formats[i] = "%8.0g"
        }
        else if (ntype[i] == "3") {
            types[i] = "int"
            formats[i] = "%8.0g"
        }
        else if (ntype[i] == "2") {
            types[i] = "long"
            formats[i] = "%12.0g"
        }
        else if (ntype[i] == "1") {
            types[i] = "float"
            formats[i] = "%9.0g"
        }
        else {
            types[i] = "double"
            formats[i] = "%10.0g"
        }
    }
    st_local("longcolumns", invtokens(longcols))
    strl = types :== "strL"
    if (mode == 1) {
        rowwidth = 0
        for (i = 1; i <= k; i++) {
            if (strl[i]) continue
            if (substr(types[i], 1, 3) == "str") rowwidth = rowwidth + strtoreal(substr(types[i], 4, .))
            else if (types[i] == "byte") rowwidth = rowwidth + 1
            else if (types[i] == "int") rowwidth = rowwidth + 2
            else if (types[i] == "double") rowwidth = rowwidth + 8
            else rowwidth = rowwidth + 4
        }
        worth = any(strl) & !anyof(names, "_all") & !anyof(names, "_skip") &
            strtoreal(st_local("strlbytes")) >= strtoreal(st_local("nobs")) * rowwidth
        if (worth & k <= 32767) {
            st_local("__cimport_fullnames", invtokens(names))
            st_local("__cimport_fulltypes", invtokens(types))
            st_local("__cimport_full", "1")
        }
        if (worth & !all(strl) & sum(strl) <= 32767 & strlen(invtokens(names)) < c("max_cmdlen") - 16) {
            st_local("__cimport_strlnames", invtokens(select(names, strl)))
            st_local("__cimport_dta", "1")
        }
        if (st_local("__cimport_full") == "1" | st_local("__cimport_dta") == "1") return
        for (i = 1; i <= k; i++) st_global("_cimport_label_" + strofreal(i), "")
    }
    if (mode == 3) {
        /* Every variable came from the file, in order. Creating variables
           marks the data changed, so a variable is added and dropped: under
           version 14 neither characteristics nor labels do, nor a variable
           with a temporary name (about a millisecond on 75,000 rows). */
        for (i = 1; i <= k; i++) {
            st_varformat(i, formats[i])
            if (labels[i] != "") st_varlabel(i, labels[i])
        }
        tag = "cimport_changed"
        for (suffix = 1; _st_varindex(tag) < .; suffix++) tag = "cimport_changed" + strofreal(suffix)
        st_dropvar(st_addvar("byte", tag))
        return
    }
    if (mode == 2) {
        /* The strL variables are in memory, in file order. */
        nostrl = !strl
        (void) st_addvar(select(types, nostrl), select(names, nostrl))
        stata("order " + invtokens(names))
        index = st_varindex(names)
    }
    else index = st_addvar(types, names)
    for (i = 1; i <= k; i++) {
        st_varformat(index[i], formats[i])
        if (labels[i] != "") st_varlabel(index[i], labels[i])
    }
}
/* Native infers display formats from a sample of the file's lines: the
   first 513 (header included), then the later lines in groups of 512 counted
   from the second line, group b when b is a power of two from 2, and a last
   group that is not full. Integers there display as %8.0g, or %12.0g in long
   range; a non-integer as float (%9.0g), or double (%10.0g) with asdouble or
   an integer beyond long range in the sample. The sample is read 256 numeric
   variables at a time. */
void _cimport_sample_formats(real scalar hasheader, real scalar asdouble)
{
    real scalar j, n, b, last, first, stop
    real colvector obs
    real rowvector numvars, idx, count, nonint, nint, lo, hi
    real matrix x, isint
    n = st_nobs()
    if (n < 1) return
    obs = 1::min((n, 513 - hasheader))
    for (b = 2; 512 * b + 2 - hasheader <= n; b = 2 * b) {
        obs = obs \ ((512 * b + 2 - hasheader)::min((n, 512 * b + 513 - hasheader)))
    }
    last = n - 2 + hasheader
    if (last >= 512 & mod(last + 1, 512)) obs = obs \ ((512 * floor(last / 512) + 2 - hasheader)::n)
    obs = uniqrows(obs)
    numvars = J(1, 0, .)
    for (j = 1; j <= st_nvar(); j++) {
        if (st_isnumvar(j)) numvars = numvars, j
    }
    for (first = 1; first <= cols(numvars); first = stop + 1) {
        stop = min((first + 255, cols(numvars)))
        idx = numvars[|first \ stop|]
        x = st_data(obs, idx)
        count = colsum(x :< .)
        isint = (x :< .) :& (x :== floor(x))
        nint = colsum(isint)
        nonint = count - nint
        x = editmissing(x, 0)
        lo = colmin(isint :* x + (!isint) :* maxdouble())
        hi = colmax(isint :* x - (!isint) :* maxdouble())
        for (j = 1; j <= cols(idx); j++) {
            if (!count[j]) continue
            if (!nint[j]) lo[j] = hi[j] = 0
            if (nonint[j]) st_varformat(idx[j], asdouble | lo[j] < -2147483647 | hi[j] > 2147483620 ? "%10.0g" : "%9.0g")
            else if (lo[j] >= -32767 & hi[j] <= 32740) st_varformat(idx[j], "%8.0g")
            else if (lo[j] >= -2147483647 & hi[j] <= 2147483620) st_varformat(idx[j], "%12.0g")
            else st_varformat(idx[j], "%10.0g")
        }
    }
}
/* 1 if Mata can open the file for reading. */
real scalar _cimport_canopen(string scalar path)
{
    real scalar fh
    fh = _fopen(path, "r")
    if (fh < 0) return(0)
    fclose(fh)
    return(1)
}
/* The next values of the strL transfer file: their lengths, then their bytes
   cut with one substr(). */
string colvector _cimport_blob_cut(real scalar fh, real colvector len)
{
    real scalar total
    string scalar buf
    total = sum(len)
    if (!total) return(J(rows(len), 1, ""))
    buf = fread(fh, total)
    if (strlen(buf) != total) _error(692, "cimport: damaged strL transfer file")
    return(substr(buf, runningsum(len) - len :+ 1, len))
}
/* Stores the strL values of cimport blobmem (src/cimport/cimport_impl.c), a
   block of rows at a time, with nothing written to a file. For each block,
   cimport blobfill stores the block's nonempty values in a frame of this
   function's: their rows, positions and lengths, and their kind (0 the
   source bytes, 1 with "" undone to ", 2 with CR LF undone to CR, 3 both, 8
   bytes of the block's text). Mata reads the block's byte range of the
   source once and cuts the values from it, or from the text, which arrives
   in global macros. A column more than half full lists every row of the
   block (empty ones with length 0) and is stored as one range; otherwise its
   empty values are not stored, as the variables were created with every
   observation empty. Frames need Stata 16; the
   frame commands are run as strings so that this file still compiles in
   Stata 14 and 15, which use _cimport_read_blobs instead. */
void _cimport_read_blobmem(string scalar path, real rowvector vars, real scalar n)
{
    real scalar src, nb, most, g, row, b, start, rlen, ntext, nent, k, c, off, at, rc, piece
    real rowvector info, counts
    real matrix X, E
    real colvector t, sel
    string scalar region, text, fr, home, q, cmd
    string colvector parts, s

    nb = strtoreal(st_global("_cimport_blobblocks"))
    most = strtoreal(st_global("_cimport_blobmax"))
    if (nb >= . | most >= .) _error(692, "cimport: damaged strL transfer")
    if (!most) return
    src = _fopen(path, "r")
    if (src < 0) _error(603, "cimport: could not reopen the file to read its strL values")
    fseek(src, 0, 1)
    if (ftell(src) != strtoreal(st_global("_cimport_blobsource"))) {
        fclose(src)
        _error(692, "cimport: the file changed while it was imported")
    }
    q = char(34)
    piece = max((1024, c("macrolen") - 64))
    home = c("frame")
    fr = st_tempname()
    rc = _stata("version 16: frame create " + fr)
    if (!rc) rc = _stata("version 16: frame " + fr + ": quietly set obs " + strofreal(most, "%20.0f"))
    if (!rc) rc = _stata("version 16: frame " + fr + ": quietly generate long r = .")
    if (!rc) rc = _stata("version 16: frame " + fr + ": quietly generate long p = .")
    if (!rc) rc = _stata("version 16: frame " + fr + ": quietly generate long n = .")
    if (!rc) rc = _stata("version 16: frame " + fr + ": quietly generate byte t = .")
    row = 1
    for (g = 0; g < nb & !rc; g++) {
        cmd = "version 16: frame " + fr + ": plugin call ctools_plugin r p n t, " + q + "cimport " +
              st_local("threads_code") + " blobfill " + strofreal(g) + q
        rc = _stata(cmd)
        if (rc) break
        info = strtoreal(tokens(st_global("_cimport_blobinfo")))
        counts = strtoreal(tokens(st_global("_cimport_blobcounts")))
        if (cols(info) != 5 | cols(counts) != cols(vars)) { rc = 692; break; }
        b = info[1]; start = info[2]; rlen = info[3]; ntext = info[4]; nent = info[5]
        if (b < 1 | row + b - 1 > n | sum(counts) != nent) { rc = 692; break; }
        if (nent) {
            /* The current frame is changed only to copy the block. */
            rc = _stata("version 16: frame change " + fr)
            if (rc) break
            X = st_data((1, nent), (1, 2, 3, 4))
            rc = _stata("version 16: frame change " + home)
            if (rc) break
        }
        region = ""
        if (rlen) {
            fseek(src, start, -1)
            region = fread(src, rlen)
            if (strlen(region) != rlen) { rc = 692; break; }
        }
        text = ""
        if (ntext) {
            parts = J(0, 1, "")
            for (at = 0; at < ntext; at = strtoreal(st_global("_cimport_blobnext"))) {
                rc = _stata("plugin call ctools_plugin, " + q + "cimport blobtext " + strofreal(g) + " " +
                            strofreal(at, "%20.0f") + " " + strofreal(piece) + q)
                if (rc) break
                parts = parts \ st_global("_cimport_blobtext")
            }
            if (rc) break
            text = invtokens(parts', "")
            if (strlen(text) != ntext) { rc = 692; break; }
        }
        off = 0
        for (k = 1; k <= cols(vars); k++) {
            c = counts[k]
            if (!c) continue
            E = X[|off + 1, 1 \ off + c, 4|]
            off = off + c
            t = E[., 4]
            s = substr(region, E[., 2], E[., 3])
            sel = selectindex(t :== 8)
            if (length(sel)) s[sel] = substr(text, E[sel, 2], E[sel, 3])
            sel = selectindex(t :== 1 :| t :== 3)
            if (length(sel)) s[sel] = subinstr(s[sel], q + q, q, .)
            sel = selectindex(t :== 2 :| t :== 3)
            if (length(sel)) s[sel] = subinstr(s[sel], char(13) + char(10), char(13), .)
            /* A column more than half full lists every row, in order. */
            if (c == b) st_sstore(row::row + b - 1, vars[k], s)
            else st_sstore((row - 1) :+ E[., 1], vars[k], s)
        }
        row = row + b
    }
    fclose(src)
    if (c("frame") != home) (void) _stata("version 16: frame change " + home)
    (void) _stata("version 16: frame drop " + fr)
    if (rc) _error(rc, "cimport: the strL transfer failed")
    if (row - 1 != n) _error(692, "cimport: damaged strL transfer")
}
/* Reads the strL transfer file of cimport_blobs (src/cimport/cimport_impl.c),
   which describes its layout. Per block of rows, the block's byte range of the
   source file is read once; each strL variable's values are cut from it with
   one substr(), apart from values decoding changed, which follow as bytes;
   then each variable's block is stored with one st_sstore(). Empty values
   are not stored when the file lists the nonempty rows: the variables were
   created with every observation empty. */
void _cimport_read_blobs(string scalar path, real rowvector vars, real scalar n)
{
    real scalar fh, src, k, row, b, start, rlen, mode, p, format, ncols, direct, size, plen, m
    real colvector io, pos, len, idx, rows
    string scalar region, source
    string colvector s
    fh = fopen(path, "r")
    io = bufio()
    format = fbufget(io, fh, "%4bu")
    ncols = fbufget(io, fh, "%4bu")
    direct = fbufget(io, fh, "%4bu")
    size = fbufget(io, fh, "%8z")
    plen = fbufget(io, fh, "%4bu")
    if (format != 3 | ncols != cols(vars)) _error(692, "cimport: damaged strL transfer file")
    src = -1
    if (direct) {
        source = plen ? fread(fh, plen) : ""
        if (strlen(source) != plen) _error(692, "cimport: damaged strL transfer file")
        src = _fopen(source, "r")
        if (src < 0) _error(603, "cimport: could not reopen the file to read its strL values")
        fseek(src, 0, 1)
        if (ftell(src) != size) _error(692, "cimport: the file changed while it was imported")
    }
    for (row = 1; row <= n; row = row + b) {
        b = fbufget(io, fh, "%4bu")
        if (b < 1 | row + b - 1 > n) _error(692, "cimport: damaged strL transfer file")
        start = fbufget(io, fh, "%8z")
        rlen = fbufget(io, fh, "%8z")
        region = ""
        if (rlen) {
            if (src < 0) _error(692, "cimport: damaged strL transfer file")
            fseek(src, start, -1)
            region = fread(src, rlen)
            if (strlen(region) != rlen) _error(692, "cimport: the file changed while it was imported")
        }
        for (k = 1; k <= cols(vars); k++) {
            mode = fbufget(io, fh, "%4bu")
            if (mode == 2) continue
            m = b
            rows = row::row + b - 1
            if (mode >= 4) {
                mode = mode - 4
                m = fbufget(io, fh, "%4bu")
                if (m < 1 | m > b) _error(692, "cimport: damaged strL transfer file")
                idx = fbufget(io, fh, "%4bu", m, 1)
                if (min(idx) < 1 | max(idx) > b) _error(692, "cimport: damaged strL transfer file")
                rows = (row - 1) :+ idx
            }
            if (mode == 0) {
                s = _cimport_blob_cut(fh, fbufget(io, fh, "%4bu", m, 1))
            }
            else if (mode == 1) {
                pos = fbufget(io, fh, "%4bu", m, 1)
                len = fbufget(io, fh, "%4bu", m, 1)
                s = substr(region, pos, len)
                p = fbufget(io, fh, "%4bu")
                if (p) {
                    idx = fbufget(io, fh, "%4bu", p, 1)
                    if (min(idx) < 1 | max(idx) > m) _error(692, "cimport: damaged strL transfer file")
                    s[idx] = _cimport_blob_cut(fh, fbufget(io, fh, "%4bu", p, 1))
                }
            }
            else _error(692, "cimport: damaged strL transfer file")
            st_sstore(rows, vars[k], s)
        }
    }
    if (src >= 0) fclose(src)
    fclose(fh)
}
end

* Native row/column ranges support f:l as well as omitted bounds.
program define _cimport_range, rclass
    version 14.1
    * kind is row (rowrange()) or column (colrange()); messages as native's
    args range kind
    local opt = cond("`kind'" == "column", "colrange", "rowrange")
    local short = cond("`kind'" == "column", "col", "row")
    local colon = strpos("`range'",":")
    local first 0
    local last 0
    if `colon' {
        local left = strtrim(substr("`range'",1,`colon'-1))
        local right = strtrim(substr("`range'",`colon'+1,.))
        if "`left'"=="f" local left "1"
        if "`left'"!="" {
            local first = real("`left'")
            if missing(`first') | `first'<1 | strpos("`left'",".") {
                di as error "{bf:`left'} is an invalid start `kind' in {bf:`opt'()} option"
                exit 198
            }
        }
        if inlist("`right'","","l") & "`left'"!="" local last 0
        else {
            local last = real("`right'")
            if missing(`last') | `last'<1 | (`colon'==1 & strpos("`right'",".")) {
                di as error "{bf:`right'}: invalid end `short' in {bf:`opt'()} option"
                exit 198
            }
        }
        if `last' & `first'>`last' {
            di as error "{bf:`range'} is an invalid `kind' range in {bf:`opt'()} option"
            di as error "start `kind' is greater than end `kind'"
            exit 198
        }
    }
    else {
        local first = real("`range'")
        if missing(`first') | `first'<1 | strpos("`range'",".") {
            di as error "{bf:`range'} is an invalid start `kind' in {bf:`opt'()} option"
            exit 198
        }
    }
    return scalar first=`first'
    return scalar last=`last'
end

* The programs below were separate ado files until ctools 1.3.2. Stata's
* net install takes at most 100 files per package, so helpers that only
* this file calls are kept here as its subprograms.

* Binary statistical readers; parsing and values are provided by C.
program define _cio_import, rclass
    version 14.1
    gettoken kind 0 : 0, parse(" ,")
    local options CLEAR SNAPshot FAST
    if inlist("`kind'","sas","spss","sasxport8","dbase","xls","xlsx") local options "`options' CASE(string)"
    if "`kind'"=="sas" local options "`options' ENCoding(string) BCAT(string)"
    if "`kind'"=="sasxport8" local options "`options' VLABFile(string)"
    if "`kind'"=="spss" local options "`options' ENCoding(string) ZSAV"
    if "`kind'"=="sasxport5" local options "`options' Describe NOVALlabels Member(string)"
    if inlist("`kind'","xls","xlsx") local options "`options' SHEET(string) CELLRange(string) FIRSTrow ALLString ALLSTRINGFMT(string) DEScribe"
    * Parse selectors as text: their variables belong to the incoming dataset.
    syntax [anything(everything equalok)] [, `options']
    if "`kind'"=="sasxport5" & "`describe'"!="" & ("`clear'"!="" | "`novallabels'"!="") exit 198
    if "`snapshot'" != "" & "`fast'" != "" {
        di as error "options snapshot and fast may not be combined"
        exit 198
    }
    local rest `"`macval(anything)'"'
    local prefix ""
    local found 0
    while `"`macval(rest)'"'!="" {
        gettoken token rest : rest, bind quotes
        if `"`macval(token)'"'=="using" {
            local found 1
            continue, break
        }
        local prefix `"`macval(prefix)' `macval(token)'"'
    }
    if `found' {
        gettoken using extra : rest
        local selection `"`macval(prefix)' `macval(extra)'"'
    }
    else {
        gettoken using selection : anything
    }
    if `"`macval(using)'"'=="" exit 198
    local rest = strtrim(`"`macval(selection)'"')
    local wanted ""
    local selectors ""
    while `"`macval(rest)'"'!="" {
        gettoken token rest : rest, bind quotes
        if inlist(`"`macval(token)'"',"if","in") {
            local selectors `"`macval(token)' `macval(rest)'"'
            continue, break
        }
        local wanted `"`macval(wanted)' `macval(token)'"'
    }
    local wanted = strtrim(`"`wanted'"')
    if !inlist("`kind'","sas","spss","xls","xlsx") & `"`wanted'"'!="" exit 198
    if !inlist("`kind'","sas","spss") & `"`selectors'"'!="" exit 198
    * SAS/SPSS: the reader applies the namelist and a plain in range while it
    * reads, as native does, so only the selection is ever held in memory.
    * Native reads the namelist as exact (case-sensitive) names, with a
    * hyphen separating names as a blank does; names the file lacks are
    * ignored. The if is applied after the read, to the selected variables.
    local __cio_select ""
    local __cio_rowfirst 0
    local __cio_rowlast 0
    if inlist("`kind'","sas","spss") {
        local __cio_select : subinstr local wanted "-" " ", all
        local rest `"`macval(selectors)'"'
        local kept ""
        while `"`macval(rest)'"' != "" {
            gettoken token rest : rest, bind quotes
            if `"`token'"' == "in" & !`__cio_rowlast' {
                gettoken range after : rest, bind
                if regexm(`"`range'"', "^([fF]|[0-9]+)(/([0-9]+))?$") {
                    local first = regexs(1)
                    local last = regexs(3)
                    if inlist("`first'","f","F") local first 1
                    if "`last'" == "" local last `first'
                    if `first' >= 1 & `last' >= `first' {
                        local __cio_rowfirst `first'
                        local __cio_rowlast `last'
                        local rest `"`macval(after)'"'
                        continue
                    }
                }
            }
            local kept `"`macval(kept)' `macval(token)'"'
        }
        local selectors = strtrim(`"`macval(kept)'"')
    }
    local ext ".sas7bdat"
    if "`kind'" == "xls" local ext ".xls"
    if "`kind'" == "xlsx" local ext ".xlsx"
    if "`kind'" == "shp" local ext ".shp"
    if "`kind'" == "dbase" local ext ".dbf"
    if "`kind'" == "spss" local ext ".sav"
    if "`kind'" == "spss" & "`zsav'" != "" local ext ".zsav"
    if "`kind'" == "sasxport5" local ext ".xpt"
    if "`kind'" == "sasxport8" local ext ".v8xpt"
    * ~ is the home directory, as for Stata's own file commands
    _ctools_home using `"`macval(using)'"'
    mata: st_local("base",pathbasename(st_local("using")))
    if !strpos(`"`macval(base)'"',".") {
        local using `"`macval(using)'`ext'"'
    }
    confirm file `"`macval(using)'"'
    if "`case'" == "" local case = cond("`kind'"=="sasxport5","lower","preserve")
    if "`case'"!="" {
        if substr("lower",1,strlen("`case'"))=="`case'" local case "lower"
        if substr("upper",1,strlen("`case'"))=="`case'" local case "upper"
        if strlen("`case'")>=3 & substr("preserve",1,strlen("`case'"))=="`case'" local case "preserve"
    }
    if inlist("`case'","low","u","pre") {
        if "`case'" == "low" local case "lower"
        if "`case'" == "u" local case "upper"
        if "`case'" == "pre" local case "preserve"
    }
    if !inlist("`case'","preserve","lower","upper") exit 198
    * error 4 shows Stata's own message ("no; data in memory would be lost")
    if "`clear'" == "" & "`describe'" == "" {
        if "`kind'"=="sasxport5" {
            if c(changed) error 4
        }
        else if (_N | c(k)) error 4
    }
    if `"`bcat'"' != "" {
        mata: st_local("base",pathbasename(st_local("bcat")))
        if !strpos(`"`base'"',".") local bcat `"`bcat'.sas7bcat"'
        confirm file `"`bcat'"'
    }
    if "`kind'"=="sasxport8" & `"`vlabfile'"'!="" {
        mata: st_local("vlabbase",pathbasename(st_local("vlabfile")))
        if !strpos(`"`vlabbase'"',".") local vlabfile `"`vlabfile'.v8xpt"'
        confirm file `"`vlabfile'"'
    }
    * Explicit snapshot protects the caller before any destructive operation.
    * preserve is program-scoped: early scan/load errors and Break restore it.
    local snap 0
    local snapframe ""
    if "`snapshot'" != "" & "`describe'" == "" {
        preserve
        local snap 1
    }
    * The plugin is checked before the data are cleared.
    _ctools_load
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    * Native statistical/DBF/SHP readers clear before parsing the file, even
    * when parsing fails. Excel validates its container before replacing data.
    if "`describe'"=="" & !inlist("`kind'","xls","xlsx") quietly clear
    local __cio_member "`member'"
    local __cio_sheet `"`macval(sheet)'"'
    local __cio_cellrange "`cellrange'"
    local __cio_firstrow = ("`firstrow'" != "")
    local __cio_maxvar = c(maxvar)
    local __cio_maxmemory = c(max_memory)
    local __cio_allstring = ("`allstring'" != "")
    local __cio_allstringfmt `"`allstringfmt'"'
    local __cio_filename `"`macval(using)'"'
    local __cio_encoding `"`encoding'"'
    local __cio_bcat `"`bcat'"'
    local __cio_xpf ""
    local __cio_rawstorage 0
    if "`kind'" == "sasxport5" & "`novallabels'" == "" & "`describe'" == "" {
        * The file's directory with its trailing separator ("" in the working
        * directory), from pathbasename(), which Stata 14.1 has.
        mata: st_local("xpfdir",substr(st_local("using"),1,strlen(st_local("using"))-strlen(pathbasename(st_local("using")))))
        foreach candidate in formats.xpf FORMATS.xpf {
            capture confirm file `"`xpfdir'`candidate'"'
            if !_rc local __cio_xpf `"`xpfdir'`candidate'"'
        }
    }
    if "`describe'"!="" & inlist("`kind'","xls","xlsx") {
        capture noisily plugin call ctools_plugin, "cio workbook `kind'"
        local rc=_rc
        if `rc' {
            macro drop CTOOLS_CIO_*
            exit `rc'
        }
        local sheets $CTOOLS_CIO_nworksheets
        return scalar N_worksheet=`sheets'
        * Sheet names come from the file: they reach r() and the screen
        * without passing through macro expansion.
        forvalues j=1/`sheets' {
            mata: st_local("worksheet",st_global("CTOOLS_CIO_worksheet_"+st_local("j")))
            mata: st_local("range",st_global("CTOOLS_CIO_range_"+st_local("j")))
            return local worksheet_`j' : copy local worksheet
            return local range_`j' : copy local range
            mata: printf("%s  %s\n", st_local("worksheet"), st_local("range"))
        }
        macro drop CTOOLS_CIO_*
        exit
    }
    if "`kind'"=="sasxport8" & `"`vlabfile'"'!="" {
        local __cio_filename `"`vlabfile'"'
        local __cio_rawstorage 1
        capture noisily plugin call ctools_plugin, "cio scan sasxport8"
        local rc=_rc
        if `rc' {
            capture plugin call ctools_plugin, "cio clear"
            macro drop CTOOLS_CIO_*
            exit `rc'
        }
        * The native command loads its companion before reading the main file.
        * Its second read fails before compression or case conversion, leaving
        * the companion's declared storage and original names in memory.
        if $CTOOLS_CIO_nvar>0 {
            local nobs $CTOOLS_CIO_nobs
            local nvar $CTOOLS_CIO_nvar
            local nlabels $CTOOLS_CIO_nlabels
            capture noisily _cio_import_load sasxport8 `nobs' `nvar' `nlabels', case(preserve)
            local rc=_rc
            capture plugin call ctools_plugin, "cio clear"
            macro drop CTOOLS_CIO_*
            if `rc' exit `rc'
            return scalar N=_N
            return scalar k=c(k)
            return scalar width=c(width)
            return scalar changed=c(changed)
            local datalabel : data label
            return local datalabel : copy local datalabel
            di as error "no; dataset in memory has changed since last saved"
            exit 4
        }
        local __cio_rawstorage 0
        local __cio_filename `"`macval(using)'"'
    }
    capture noisily plugin call ctools_plugin, "cio scan `kind'"
    local rc = _rc
    if `rc' {
        macro drop CTOOLS_CIO_*
        exit `rc'
    }
    local nobs = $CTOOLS_CIO_nobs
    local nvar = $CTOOLS_CIO_nvar
    local nlabels = $CTOOLS_CIO_nlabels
    if "`describe'" != "" {
        local nmembers $CTOOLS_CIO_nmembers
        * Member names, variable names and labels come from the file: they are
        * handled in Mata and reach r() and the screen without macro expansion.
        mata: st_local("members", strlower(st_global("CTOOLS_CIO_members")))
        mata: st_local("descriptions", st_local("member"))
        if "`kind'"=="sasxport5" & "`member'"=="" local descriptions : copy local members
        mata: st_local("ndescriptions", strofreal(cols(tokens(st_local("descriptions")))))
        * An error or Break in the listing still releases the scanned file
        capture noisily {
            forvalues m=1/`ndescriptions' {
                mata: st_local("item", tokens(st_local("descriptions"))[`m'])
                if "`kind'"=="sasxport5" {
                    local __cio_member : copy local item
                    plugin call ctools_plugin, "cio scan `kind'"
                    local nobs $CTOOLS_CIO_nobs
                    local nvar $CTOOLS_CIO_nvar
                    local size $CTOOLS_CIO_size
                }
                mata: printf("%s: %s variables, %s observations\n", st_local("item"), st_local("nvar"), st_local("nobs"))
                forvalues j=1/`nvar' {
                    plugin call ctools_plugin, "cio column `j'"
                    mata: printf("%s  %s  %s\n", st_global("CTOOLS_CIO_name"), st_global("CTOOLS_CIO_type"), st_global("CTOOLS_CIO_label"))
                }
            }
        }
        local rc = _rc
        capture plugin call ctools_plugin, "cio clear"
        macro drop CTOOLS_CIO_*
        if `rc' exit `rc'
        return scalar N = `nobs'
        return scalar k = `nvar'
        if "`kind'"=="sasxport5" {
            return scalar size = `size'
            if "`member'"=="" {
                return scalar n_members = `nmembers'
                return local members : copy local members
            }
        }
        exit
    }
    * A failure after the scan (allocation, SPI, Break) leaves the data in
    * memory as they were. Only an Excel workbook reaches this point with data
    * in memory (the other readers cleared them above, as native ones do);
    * from Stata 16 they are copied into a frame of this program's own, in
    * memory (the copy fails before anything changes if memory is short), and
    * restored from it with frame copy (see cimport.ado for why the load is
    * not built in a frame that is then renamed). Option snapshot preserves
    * the data instead (on disk above c(max_preservemem)), as do loads before
    * Stata 16 and loads into an empty dataset. Option fast keeps no copy and
    * undoes nothing. preserve and the frame belong to this program, so these
    * steps are inline.
    if "`fast'" == "" & !`snap' {
        if "`snapshot'" == "" & c(stata_version) >= 16 & (_N > 0 | c(k) > 0) {
            tempname snapframe
            capture noisily frame copy `c(frame)' `snapframe'
        }
        else {
            capture noisily preserve
            local snap 1
        }
        * A failed copy still releases the scanned file
        if _rc {
            local rc = _rc
            capture plugin call ctools_plugin, "cio clear"
            macro drop CTOOLS_CIO_*
            exit `rc'
        }
    }
    local loadcase "`case'"
    if inlist("`kind'","sas","spss") | (inlist("`kind'","xls","xlsx") & "`firstrow'"=="") local loadcase preserve
    capture noisily _cio_import_load `kind' `nobs' `nvar' `nlabels', case(`loadcase') wanted(`"`wanted'"') allstringfmt(`"`allstringfmt'"') `novallabels' `allstring' filename(`"`macval(using)'"')
    local rc = _rc
    capture plugin call ctools_plugin, "cio clear"
    macro drop CTOOLS_CIO_*
    if !`rc' {
        capture noisily {
            if `"`selectors'"' != "" quietly keep `selectors'
            if inlist("`kind'","sas","spss") {
                if "`case'" != "preserve" rename _all, `case'
                quietly compress
            }
            * Validation hook: fail as Break would once the data are loaded
            if "$CTOOLS_cimport_fault" == "load" error 1
        }
        local rc = _rc
    }
    if `rc' {
        nobreak {
            if "`snapframe'" != "" frame copy `snapframe' `c(frame)', replace
            else if `snap' restore
        }
        exit `rc'
    }
    if "`snapframe'" != "" frame drop `snapframe'
    else if `snap' restore, not
    return scalar N = _N
    return scalar k = c(k)
    return scalar width = c(width)
    return scalar changed = 1
    local datalabel : data label
    return local datalabel : copy local datalabel
    return local filename `"`macval(using)'"'
end

program define _cio_import_load
    version 14.1
    gettoken kind 0 : 0, parse(" ,")
    gettoken nobs 0 : 0, parse(" ,")
    gettoken nvar 0 : 0, parse(" ,")
    gettoken nlabels 0 : 0, parse(" ,")
    syntax, CASE(string) [WANTED(string) NOVALLabels ALLSTRINGFMT(string) ALLString FILEname(string)]
    * Plugin registrations belong to the invoking ado program.
    _ctools_load
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    clear
    * The file's data label. label data cannot quote one holding an
    * unbalanced compound quote (it fails or stores other text); that label
    * comes from an empty dataset written to hold it. Later steps keep it.
    mata: st_local("__data_label",st_global("CTOOLS_CIO_datalabel"))
    capture label data `"`macval(__data_label)'"'
    local __label_rc = _rc
    local __label_got : data label
    mata: st_local("__label_rc",strofreal(`__label_rc' | st_local("__label_got")!=usubstr(st_local("__data_label"),1,80)))
    if `__label_rc' {
        tempfile __label_stub
        mata: _cio_label_stub(st_local("__label_stub"),st_local("__data_label"))
        quietly use `"`__label_stub'"', clear
    }
    mata: __ctools_cio_sets = J(1,0,"")
    if `nvar' {
    if inlist("`kind'","xls","xlsx") {
        * Resolve the entire header before creating variables: a later column's
        * fallback can require moving an earlier header to its own column letter.
        forvalues j=1/`nvar' {
            plugin call ctools_plugin, "cio column `j'"
            mata: st_local("__cio_candidate_`j'",_cio_excel_name(st_global("CTOOLS_CIO_name"),st_global("CTOOLS_CIO_excel_column"),st_local("case")))
            mata: st_local("__cio_fallback_`j'",st_global("CTOOLS_CIO_excel_column"))
        }
        mata: _cio_excel_names(`nvar')
    }
    forvalues j=1/`nvar' {
        plugin call ctools_plugin, "cio column `j'"
        if inlist("`kind'","xls","xlsx") {
            local vname "`__cio_resolved_`j''"
        }
        else mata: st_local("vname",ustrtoname(st_global("CTOOLS_CIO_name")))
        if "`case'" == "lower" & !inlist("`kind'","xls","xlsx") local vname = ustrlower("`vname'")
        if "`case'" == "upper" & !inlist("`kind'","xls","xlsx") local vname = ustrupper("`vname'")
        if inlist("`kind'","xls","xlsx") {
            if "`case'" == "lower" local vname = ustrlower("`vname'")
            if "`case'" == "upper" local vname = ustrupper("`vname'")
        }
        if !inlist("`kind'","xls","xlsx") {
            capture confirm new variable `vname'
            if _rc local vname "v`j'"
        }
        if $CTOOLS_CIO_external local external `external' `j'
        if $CTOOLS_CIO_numerictext local textnumbers `textnumbers' `j'
        local sourcecol`j' "${CTOOLS_CIO_excel_column}"
        local sourcevar`j' "`vname'"
        local type ${CTOOLS_CIO_type}
        if inlist("`kind'","xls","xlsx") & ("`vname'"=="_cons" | ustrregexm("`vname'","^str[0-9]+$")) {
            di as error "`vname' invalid varname"
            exit 198
        }
        mata: (void) st_addvar(st_local("type"),st_local("vname"))
        mata: st_varlabel(st_local("vname"),st_global("CTOOLS_CIO_labelnumber")=="" ? st_global("CTOOLS_CIO_label") : strtrim(strofreal(strtoreal(st_global("CTOOLS_CIO_labelnumber")),"%18.0g")))
        mata: st_varformat(st_local("vname"),st_global("CTOOLS_CIO_format"))
        if ("`novallabels'" == "" | "`kind'"=="sasxport5") & inlist("`kind'","spss","sasxport5") {
            * Label-set names come from the file and stay in Mata. Native
            * SPSS labels use their dictionary ordinal, beginning at 0.
            mata: st_local("labelname",_cio_labelset_name(st_local("kind"),st_global("CTOOLS_CIO_labelset"),1))
            * String variables take no value label: native SPSS notes it,
            * native XPORT fails unless novallabels is given.
            if "`labelname'" != "" & substr("`type'",1,3) == "str" {
                if "`kind'" == "sasxport5" & "`novallabels'" == "" {
                    mata: printf("{err}%s has string variable %s with format (value label) %s; value labels for string variables are not allowed in Stata\n", st_local("filename"), st_local("vname"), st_global("CTOOLS_CIO_labelset"))
                    exit 610
                }
                if "`kind'" == "spss" di as text "note: may not attach value label to string variable: no label applied"
                local labelname ""
            }
            if "`labelname'" != "" mata: st_varvaluelabel(st_local("vname"),st_local("labelname"))
        }
    }
    }
    * The observations are added after the variables: a variable added to
    * `nobs' observations is filled, and the dataset rewritten, once per
    * variable; set obs fills them all in one pass.
    quietly set obs `=max(1,`nobs')'
    if `nobs' {
    foreach hdridx of local external {
        forvalues i=1/`nobs' {
            local offset 0
            mata: __ctools_cio_blob = ""
            while 1 {
                plugin call ctools_plugin, "cio blob `hdridx' `i' `offset'"
                mata: __ctools_cio_blob = __ctools_cio_blob + _cio_hex_bytes(st_global("CTOOLS_CIO_blobhex"))
                local offset = $CTOOLS_CIO_blobnext
                if !`offset' continue, break
            }
            mata: st_sstore(`i',`hdridx',__ctools_cio_blob)
        }
        mata: mata drop __ctools_cio_blob
    }
    }
    if "`kind'" == "shp" {
        foreach field in type xmin ymin xmax ymax zmin zmax mmin mmax num_records {
            mata: st_global("_dta[shp_"+st_local("field")+"]",st_global("CTOOLS_CIO_shp_"+st_local("field")))
        }
        char _dta[dta_shp_ver] "1"
        char _dta[shp_code] "9994"
        char _dta[shp_version] "1000"
    }
    if "`kind'"=="dbase" {
        foreach field in version year month day {
            mata: st_global("_dta[dbf_"+st_local("field")+"]",st_global("CTOOLS_CIO_dbf_"+st_local("field")))
        }
    }
    if `nobs' & `nvar' plugin call ctools_plugin *, "cio load"
    * Numbers in string columns: allstring()'s format, else native's %18.0g
    * (for which the plugin reserved 18 characters).
    local textfmt `"`allstringfmt'"'
    if `"`textfmt'"' == "" local textfmt "%18.0g"
    if `nobs' {
        foreach j of local textnumbers {
            tempvar rawnumber
            quietly gen double `rawnumber' = .
            plugin call ctools_plugin `rawnumber', "cio textnumbers `j'"
            local target `sourcevar`j''
            mata: _cio_format_text(st_local("target"),st_local("rawnumber"),st_local("textfmt"),st_local("allstringfmt")=="")
            quietly drop `rawnumber'
            quietly compress `target'
            mata: st_varformat(st_local("target"),"%"+strofreal(st_vartype(st_local("target"))=="strL" ? 9 : max((9,max(strlen(st_sdata(.,st_local("target")))))))+"s")
        }
    }
    if !`nobs' quietly drop in 1
    if "`kind'" == "shp" {
        capture confirm variable shape_order
        if !_rc sort _ID
    }
    if "`novallabels'" == "" & `nlabels' {
        forvalues i=1/`nlabels' {
            plugin call ctools_plugin, "cio label `i'"
            mata: st_local("labelname",_cio_labelset_name(st_local("kind"),st_global("CTOOLS_CIO_labelset"),0))
            mata: st_vlmodify(st_local("labelname"),strtoreal(st_global("CTOOLS_CIO_labelvalue")),st_global("CTOOLS_CIO_labeltext"))
        }
    }
    if inlist("`kind'","xls","xlsx") & `"`wanted'"'!="" {
        local rest `"`wanted'"'
        local desired ""
        local selected ""
        local position 0
        while strtrim(`"`rest'"')!="" {
            gettoken name rest : rest, parse(" =")
            if `: list name in desired' exit 198
            local ++position
            gettoken equal tail : rest, parse(" =")
            local column ""
            if "`equal'"=="=" {
                gettoken column rest : tail
                local column=upper("`column'")
            }
            else {
                mata: st_local("column",_cio_excel_column(`position'))
            }
            local index 0
            forvalues j=1/`nvar' {
                if "`sourcecol`j''"=="`column'" local index `j'
            }
            tempvar picked`position'
            if !`index' {
                if !`nvar' exit 198
                mata: st_local("columnindex",strofreal(_cio_excel_index(st_local("column"))))
                mata: st_local("columnlimit",strofreal(_cio_excel_index(st_local("sourcecol`nvar'"))+1))
                if `columnindex'<1 | `columnindex'>`columnlimit' exit 198
                if "`allstring'"!="" quietly gen str1 `picked`position''=""
                else {
                    quietly gen byte `picked`position''=.
                    format `picked`position'' %10.0g
                }
            }
            else clonevar `picked`position''=`sourcevar`index''
            local selected `selected' `picked`position''
            local desired `desired' `name'
        }
        keep `selected'
        local j=0
        foreach v of local selected {
            local ++j
            local name : word `j' of `desired'
            rename `v' `name'
        }
    }
    capture mata: mata drop __ctools_cio_sets
end

mata:
real scalar _cio_excel_index(string scalar column) {
    real scalar j, n, code
    n=0
    for(j=1;j<=strlen(column);j++) {
        code=ascii(substr(column,j,1))-64
        if(code<1 || code>26)return(.)
        n=n*26+code
        if(n>16384)return(.)
    }
    return(n)
}
string scalar _cio_excel_column(real scalar index) {
    string scalar column
    column=""
    while(index>0) {
        column=char(65+mod(index-1,26))+column
        index=floor((index-1)/26)
    }
    return(column)
}
/* Lowercase hex pairs to bytes, decoded as whole vectors. */
string scalar _cio_hex_bytes(string scalar hex) {
    real matrix digits
    if (strlen(hex) < 2) return("")
    digits = ascii(hex)
    digits = colshape(digits :- 48 :- 39 :* (digits :> 57), 2)
    return(char((digits[., 1] * 16 + digits[., 2])'))
}
end

mata:
/* An empty .dta (format 118, no variables or observations) whose header
   holds label, the first 80 characters as label data keeps them: using it
   sets a data label that label data cannot quote (an unbalanced compound
   quote). */
void _cio_label_stub(string scalar path, string scalar label) {
    string scalar head, body, map, text
    real colvector at
    real scalar fh, i
    string rowvector tags
    text=usubstr(label,1,80)
    head="<stata_dta><header><release>118</release><byteorder>LSF</byteorder><K>"+_cio_le(0,2)+"</K><N>"+_cio_le(0,8)+"</N><label>"+_cio_le(strlen(text),2)+text+"</label><timestamp>"+char(0)+"</timestamp></header>"
    tags=("variable_types","varnames","sortlist","formats","value_label_names","variable_labels","characteristics","data","strls","value_labels")
    at=J(14,1,0)
    at[2]=strlen(head)
    body=""
    for (i=1;i<=10;i++) {
        at[i+2]=at[2]+strlen("<map>")+14*8+strlen("</map>")+strlen(body)
        body=body+"<"+tags[i]+">"+(tags[i]=="sortlist" ? _cio_le(0,2) : "")+"</"+tags[i]+">"
    }
    at[13]=at[2]+strlen("<map>")+14*8+strlen("</map>")+strlen(body)
    at[14]=at[13]+strlen("</stata_dta>")
    map="<map>"
    for (i=1;i<=14;i++) map=map+_cio_le(at[i],8)
    map=map+"</map>"
    if (fileexists(path)) unlink(path)
    fh=fopen(path,"w")
    fwrite(fh,head+map+body+"</stata_dta>")
    fclose(fh)
}
string scalar _cio_le(real scalar value, real scalar bytes) {
    string scalar out
    real scalar i, rest
    out=""
    rest=value
    for (i=1;i<=bytes;i++) {
        out=out+char(mod(rest,256))
        rest=floor(rest/256)
    }
    return(out)
}
void _cio_format_text(string scalar target, string scalar source, string scalar fmt, real scalar trim) {
    real colvector values, selected
    string colvector text
    real scalar width
    values=st_data(.,source)
    selected=selectindex(values:<.)
    if (!length(selected)) return
    text=trim ? strtrim(strofreal(values[selected],fmt)) : strofreal(values[selected],fmt)
    /* Widen the variable only as far as the text needs (strL past 2045) */
    width=max(strlen(text))
    if (st_vartype(target)!="strL" & width>strtoreal(substr(st_vartype(target),4,.))) stata("quietly recast "+(width>2045 ? "strL" : "str"+strofreal(width))+" "+target)
    st_sstore(selected,target,text)
}
end

mata:
/* Native cuts the header to 32 characters, keeps Unicode letters, decimal
   digits and underscores, drops leading digits, and uses the column letter
   when the result, in the requested case, is not a valid name. */
string scalar _cio_excel_name(string scalar source, string scalar fallback, string scalar casing) {
    string scalar name, cased
    name=ustrregexra(usubstr(source,1,32),"[^\p{L}\p{Nd}_]","")
    name=ustrregexrf(name,"^\p{Nd}+","")
    cased=casing=="lower" ? ustrlower(name) : (casing=="upper" ? ustrupper(name) : name)
    /* Native does not replace _cons or str#: creating them fails. */
    if (cased=="_cons" || ustrregexm(cased,"^str[0-9]+$")) return(name)
    if (name=="" || !st_isvarname(cased)) return(fallback)
    return(name)
}
/* Value-label name for a label set named in the file: for SPSS the reader's
   labels# (the set's dictionary ordinal, as native), the lowercased name for
   the other formats. "" when a variable has no label set (add=1). */
string scalar _cio_labelset_name(string scalar kind, string scalar set, real scalar add) {
    if (add & set=="") return("")
    if (kind!="spss") return(ustrlower(ustrtoname(set)))
    return(set)
}
void _cio_excel_names(real scalar n) {
    string rowvector names, fallback
    string scalar name
    real scalar j, at, previous
    transmorphic scalar used
    names=J(1,n,""); fallback=J(1,n,"")
    used=asarray_create("string",1)
    for(j=1;j<=n;j++) fallback[j]=st_local(sprintf("__cio_fallback_%g",j))
    for(j=1;j<=n;j++) {
        name=st_local(sprintf("__cio_candidate_%g",j))
        if(asarray_contains(used,name)) name=fallback[j]
        at=j
        while(asarray_contains(used,name)) {
            previous=asarray(used,name)
            names[at]=name
            asarray(used,name,at)
            at=previous
            name=fallback[at]
        }
        names[at]=name
        asarray(used,name,at)
    }
    for(j=1;j<=n;j++) st_local(sprintf("__cio_resolved_%g",j),names[j])
}
end

* cimport parquet: a native C Parquet reader with import parquet's syntax,
* variables, messages and stored results.
program define _cimport_parquet, rclass
    version 14.1
    * As import parquet: a bare filename before the comma is the using file
    gettoken filename rest : 0, parse(" ,")
    gettoken comma : rest, parse(" ,")
    if `"`filename'"' != "" & (trim(`"`comma'"') == "," | trim(`"`comma'"') == "") {
        local 0 `"using `0'"'
    }
    local __desc 0
    capture syntax using/, DESCribe
    if !_rc local __desc 1
    else {
        syntax [anything(name=vlist id="vlist")] using/ [, clear ROWRange(string) ///
            favormemory SNAPshot FAST]
        if "`snapshot'" != "" & "`fast'" != "" {
            di as error "options snapshot and fast may not be combined"
            exit 198
        }
        if "`clear'" == "" & c(changed) error 4
    }
    _ctools_load
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    * Inputs and results pass through CTOOLS_PQ_* globals, dropped at every exit
    macro drop CTOOLS_PQ_*
    local __pq_home : environment HOME
    if c(os) == "Windows" local __pq_home : environment USERPROFILE
    capture noisily mata: _cpq_path(st_local("using"))
    if _rc {
        local rc = _rc
        macro drop CTOOLS_PQ_*
        exit `rc'
    }
    if `__desc' {
        * Stata's describe decodes the data too; unreadable data show -1 rows
        capture noisily {
            _cimport_parquet_call describe
            mata: _cpq_describe()
        }
        local rc = _rc
        local k $CTOOLS_PQ_k
        local n $CTOOLS_PQ_N
        macro drop CTOOLS_PQ_*
        if `rc' exit `rc'
        return scalar n_columns = `k'
        return scalar N_rows = `n'
        exit
    }
    * Column selection, row range, supported types and variable names, in
    * the order import parquet checks them
    capture noisily {
        _cimport_parquet_call schema
        mata: _cpq_prepare(st_local("vlist"), st_local("rowrange"))
    }
    if _rc {
        local rc = _rc
        macro drop CTOOLS_PQ_*
        exit `rc'
    }
    local nout $CTOOLS_PQ_nout
    * The undo copy, as in _cio_import: a frame from Stata 16 when there are
    * data, otherwise preserve; snapshot always preserves; fast keeps none.
    local snap 0
    local snapframe ""
    if "`fast'" == "" {
        if "`snapshot'" == "" & c(stata_version) >= 16 & (_N > 0 | c(k) > 0) {
            tempname snapframe
            capture noisily frame copy `c(frame)' `snapframe'
        }
        else {
            capture noisily preserve
            local snap 1
        }
        if _rc {
            local rc = _rc
            macro drop CTOOLS_PQ_*
            exit `rc'
        }
    }
    capture noisily {
        local nobs = $CTOOLS_PQ_last - $CTOOLS_PQ_first + 1
        if `nobs' > 2147483647 {
            di as error "ctools does not support datasets exceeding 2^31 (2.147 billion) observations"
            di as error "This is a limitation of Stata's plugin API"
            exit 920
        }
        * As import parquet: an empty file, or an end row that is not an
        * integer, fails when the rows are loaded
        if $CTOOLS_PQ_badrows | `nobs' < 1 {
            di as error "out of memory"
            exit 3900
        }
        plugin call ctools_plugin, "cparquet scan"
        quietly drop _all
        * The variables are created before the observations: a variable
        * added to `nobs' observations is filled, and the dataset rewritten,
        * once per variable; set obs fills them all in one pass.
        forvalues j = 1/`nout' {
            local type ${CTOOLS_PQ_vtype_`j'}
            local name ${CTOOLS_PQ_vname_`j'}
            * import parquet's own command, so that invalid names fail alike
            if substr("`type'", 1, 3) == "str" capture noisily quietly gen `type' `name' = ""
            else capture noisily quietly gen `type' `name' = .
            if _rc exit 3598
            if "${CTOOLS_PQ_fmt_`j'}" != "" format ${CTOOLS_PQ_fmt_`j'} `name'
        }
        mata: _cpq_labels(`nout')
        quietly set obs `nobs'
        plugin call ctools_plugin *, "cparquet load"
        forvalues j = 1/`nout' {
            if ${CTOOLS_PQ_blob_`j'} {
                tempfile blobs
                global CTOOLS_PQ_blobpath `"`blobs'"'
                plugin call ctools_plugin, "cparquet blobfile `j'"
                mata: _cpq_read_blobs(st_global("CTOOLS_PQ_blobpath"), `j', `nobs')
                capture erase `"`blobs'"'
            }
        }
        * Dictionary (categorical) columns: value labels from their entries,
        * modifying a label of the same name as import parquet does
        forvalues j = 1/`nout' {
            if ${CTOOLS_PQ_dict_`j'} {
                tempfile entries
                global CTOOLS_PQ_blobpath `"`entries'"'
                plugin call ctools_plugin, "cparquet dictfile `j'"
                mata: _cpq_dict_labels(st_global("CTOOLS_PQ_blobpath"), st_global("CTOOLS_PQ_vname_`j'"), ${CTOOLS_PQ_dict_`j'})
                label values ${CTOOLS_PQ_vname_`j'} ${CTOOLS_PQ_vname_`j'}
                capture erase `"`entries'"'
            }
        }
        * Validation hook: fail as Break would once the data are loaded
        if "$CTOOLS_cimport_fault" == "load" error 1
    }
    local rc = _rc
    capture plugin call ctools_plugin, "cparquet clear"
    macro drop CTOOLS_PQ_*
    if `rc' {
        nobreak {
            if "`snapframe'" != "" frame copy `snapframe' `c(frame)', replace
            else if `snap' restore
        }
        exit `rc'
    }
    if "`snapframe'" != "" frame drop `snapframe'
    else if `snap' restore, not
    mata: printf("{txt}(%s, %s obs)\n", `nout' == 1 ? "1 var" : strofreal(`nout', "%16.0gc") + " vars", strofreal(st_nobs(), "%16.0gc"))
    return scalar k = `nout'
    return scalar N = _N
end

* Reads the footer; a file that cannot be opened fails as in import parquet.
program define _cimport_parquet_call
    version 14.1
    args op
    _ctools_load
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    capture noisily plugin call ctools_plugin, "cparquet `op'"
    local rc = _rc
    if `rc' == 603 di as error `"file {bf:${CTOOLS_PQ_hold}} could not be opened"'
    if `rc' exit `rc'
end

mata:
/* The file as import parquet finds it: no URLs, ~ for the home directory,
   relative to the working directory, .parquet added when there is no
   extension. CTOOLS_PQ_file is the absolute path. */
void _cpq_path(string scalar f)
{
    string scalar hold, home
    hold = f
    if (pathisurl(f)) {
        errprintf("invalid filename; filename may not be specified as a URL\n")
        exit(603)
    }
    if (substr(f, 1, 1) == "~" & (strlen(f) == 1 | anyof(("/", "\"), substr(f, 2, 1)))) {
        home = st_local("__pq_home")
        if (home != "") f = home + substr(f, 2, .)
    }
    if (!pathisabs(f)) f = _cpq_resolve(pwd(), f)
    if (fileexists(f) != 1) {
        if (pathsuffix(f) == "") f = f + ".parquet"
        if (fileexists(f) != 1) {
            errprintf("file {bf:%s} not found\n", hold)
            exit(601)
        }
    }
    st_global("CTOOLS_PQ_file", f)
    st_global("CTOOLS_PQ_hold", hold)
}

/* An absolute path from the working directory and a relative path, made as
   Stata 16's path resolution makes it (which Stata 14.1 lacks): forward
   slashes, no empty, . or .. components, no trailing slash. */
string scalar _cpq_resolve(string scalar base, string scalar f)
{
    string scalar s, prefix, part
    string rowvector stack
    real scalar at
    s = subinstr(base + "/" + f, "\", "/", .)
    prefix = ""
    if (substr(s, 1, 2) == "//") prefix = "//"
    else if (substr(s, 1, 1) == "/") prefix = "/"
    else if (regexm(s, "^[A-Za-z]:")) {
        prefix = substr(s, 1, 2) + "/"
        s = substr(s, 3, .)
    }
    stack = J(1, 0, "")
    while (s != "") {
        at = strpos(s, "/")
        if (at) {
            part = substr(s, 1, at - 1)
            s = substr(s, at + 1, .)
        }
        else {
            part = s
            s = ""
        }
        if (part == "" | part == ".") continue
        if (part == "..") {
            if (cols(stack)) stack = cols(stack) > 1 ? stack[1..cols(stack) - 1] : J(1, 0, "")
            continue
        }
        stack = (stack, part)
    }
    return(prefix + invtokens(stack, "/"))
}

void _cpq_describe()
{
    real scalar i, k
    k = strtoreal(st_global("CTOOLS_PQ_k"))
    printf("{txt}Contains data from {res}%s\n", st_global("CTOOLS_PQ_file"))
    printf("{txt} Observations:        {res}%s\n", strofreal(strtoreal(st_global("CTOOLS_PQ_N")), "%16.0gc"))
    printf("{txt}    Variables:        {res}%s\n", strofreal(k, "%16.0gc"))
    printf("{txt}{hline}\n")
    printf("Column                           Type\n")
    printf("{hline}\n")
    for (i = 1; i <= k; i++) {
        printf("{p 0 48}{res}{bind:%-32s }{txt}{res}{bind:%-15s}{p_end}", abbrev(st_global("CTOOLS_PQ_name_" + strofreal(i)), 32), st_global("CTOOLS_PQ_type_" + strofreal(i)))
    }
    printf("{txt}{hline}\n")
}

real scalar _cpq_range_error(string scalar why)
{
    errprintf("invalid range in option {bf:rowrange()}")
    errprintf("{p 4 4 2}%s{p_end}", why)
    return(198)
}

/* rowrange([start][:end]) in the forms import parquet accepts, with its
   messages: f and l name the first and last rows, and a negative number
   counts from the end. Returns (start, end), end missing for the last row. */
real rowvector _cpq_rowrange(string scalar r)
{
    real scalar colon, a, b
    string scalar left, right
    a = 1
    b = .
    colon = strpos(r, ":")
    if (colon > 1) {
        left = bsubstr(r, 1, colon - 1)
        a = anyof(("f", "F"), strtrim(left)) ? 1 : strtoreal(left)
        if (a == . | strpos(left, ".") > 0) exit(_cpq_range_error("The option {bf:rowrange()} requires an integer for the starting row."))
        right = bsubstr(r, colon + 1, .)
        if (!anyof(("l", "L", ""), strtrim(right))) {
            b = strtoreal(right)
            if (b == .) exit(_cpq_range_error("The option {bf:rowrange()} requires an integer for the ending row."))
        }
    }
    else if (colon == 1) {
        right = bsubstr(r, 2, .)
        b = strtoreal(right)
        if (b == . | strpos(right, ".") > 0) exit(_cpq_range_error("The option {bf:rowrange()} requires an integer for the ending row."))
    }
    else {
        a = anyof(("f", "F"), strtrim(r)) ? 1 : strtoreal(r)
        if (a == . | strpos(r, ".") > 0) exit(_cpq_range_error("The option {bf:rowrange()} requires an integer for the starting row."))
    }
    return((a, b))
}

void _cpq_prepare(string scalar list, string scalar rowrange)
{
    real scalar k, nrows, i, j, first, last, found
    string colvector names, types, invalid
    string rowvector wanted, vnames, orig
    real rowvector sel, keep, range
    k = strtoreal(st_global("CTOOLS_PQ_k"))
    nrows = strtoreal(st_global("CTOOLS_PQ_N"))
    names = J(k, 1, "")
    types = J(k, 1, "")
    for (i = 1; i <= k; i++) {
        names[i] = st_global("CTOOLS_PQ_name_" + strofreal(i))
        types[i] = st_global("CTOOLS_PQ_type_" + strofreal(i))
    }
    /* Selection by exact name, the first column of that name; a name given
       twice is read once */
    sel = J(1, 0, .)
    if (strtrim(list) != "") {
        wanted = tokens(list)
        invalid = J(0, 1, "")
        for (i = 1; i <= cols(wanted); i++) {
            found = 0
            for (j = 1; j <= k; j++) {
                if (wanted[i] == names[j]) {
                    if (!anyof(sel, j)) sel = (sel, j)
                    found = 1
                    break
                }
            }
            if (!found) invalid = (invalid \ wanted[i])
        }
        if (cols(sel) == 0) {
            errprintf("invalid column names\n")
            exit(198)
        }
        for (i = 1; i <= rows(invalid); i++) printf("{txt}note: column {bf:%s} ignored because of invalid column name.\n", invalid[i])
    }
    else {
        if (k == 0) {
            errprintf("file contains 0 column\n")
            exit(198)
        }
        sel = 1..k
    }
    first = 1
    last = nrows
    if (strtrim(rowrange) != "") {
        range = _cpq_rowrange(rowrange)
        first = range[1]
        last = range[2]
        if (first < 0) first = nrows + first + 1
        if (last == .) last = nrows
        if (last < 0) last = nrows + last + 1
        if (first < 1 | first > nrows | last < 1 | last > nrows) exit(_cpq_range_error("The specified integers are outside of the valid range."))
        if (first > last) exit(_cpq_range_error("The ending row must be greater than the starting row."))
    }
    /* Types Stata cannot hold are skipped with a note */
    keep = J(1, 0, .)
    invalid = J(0, 1, "")
    for (i = 1; i <= cols(sel); i++) {
        if (anyof(("null", "struct", "list", "map", "duration"), types[sel[i]])) invalid = (invalid \ names[sel[i]])
        else keep = (keep, sel[i])
    }
    if (cols(keep) == 0) {
        errprintf("all columns are invalid due to unsupported data types\n")
        exit(198)
    }
    for (i = 1; i <= rows(invalid); i++) printf("{txt}note: column {bf:%s} ignored because of unsupported Parquet type.\n", invalid[i])
    orig = names[keep]'
    vnames = _cpq_varnames(orig)
    st_global("CTOOLS_PQ_cols", invtokens(strofreal(keep :- 1)))
    st_global("CTOOLS_PQ_first", strofreal(first, "%21.0g"))
    st_global("CTOOLS_PQ_last", strofreal(floor(last), "%21.0g"))
    st_global("CTOOLS_PQ_badrows", strofreal(last != floor(last)))
    st_global("CTOOLS_PQ_nout", strofreal(cols(keep)))
    for (i = 1; i <= cols(keep); i++) {
        st_global("CTOOLS_PQ_vname_" + strofreal(i), vnames[i])
        st_global("CTOOLS_PQ_orig_" + strofreal(i), orig[i])
    }
}

/* Variable names as import parquet makes them: letters, decimal digits and
   underscores are kept, leading digits are dropped, and the result is cut
   to 32 characters; an empty name becomes v#, # the column's 0-based
   position among the imported ones; and if any two names coincide, every
   column is named that way. */
string rowvector _cpq_varnames(string rowvector names)
{
    real scalar i, j
    string rowvector out
    out = J(1, cols(names), "")
    for (i = 1; i <= cols(names); i++) {
        out[i] = ustrregexra(names[i], "[^\p{L}\p{Nd}_]", "")
        out[i] = ustrregexrf(out[i], "^\p{Nd}+", "")
        out[i] = usubstr(out[i], 1, 32)
        if (out[i] == "") out[i] = "v" + strofreal(i - 1)
    }
    for (i = 1; i <= cols(out); i++) {
        for (j = i + 1; j <= cols(out); j++) {
            if (out[i] == out[j]) return("v" :+ strofreal(0..cols(out) - 1))
        }
    }
    return(out)
}

/* A variable whose name differs from its column is labeled with the column */
void _cpq_labels(real scalar n)
{
    real scalar i
    string scalar name, orig
    for (i = 1; i <= n; i++) {
        name = st_global("CTOOLS_PQ_vname_" + strofreal(i))
        orig = st_global("CTOOLS_PQ_orig_" + strofreal(i))
        if (name != orig) st_varlabel(name, usubstr(orig, 1, 80))
    }
}

void _cpq_dict_labels(string scalar path, string scalar name, real scalar n)
{
    real scalar fh, i, length
    real colvector io
    string colvector text
    fh = fopen(path, "r")
    io = bufio()
    text = J(n, 1, "")
    for (i = 1; i <= n; i++) {
        length = fbufget(io, fh, "%4bu")
        if (length) text[i] = fread(fh, length)
    }
    fclose(fh)
    st_vlmodify(name, (0::n - 1), text)
}

void _cpq_read_blobs(string scalar path, real scalar col, real scalar n)
{
    real scalar fh, i, length
    real colvector io
    fh = fopen(path, "r")
    io = bufio()
    for (i = 1; i <= n; i++) {
        length = fbufget(io, fh, "%4bu")
        st_sstore(i, col, length ? fread(fh, length) : "")
    }
    fclose(fh)
}
end
