*! version 1.2.4 20261103 github.com/mdroste/stata-ctools

program define cimport, rclass
    version 14.1

    * Parse the subcommand
    gettoken subcmd 0 : 0, parse(" ,")
    if inlist("`subcmd'", "delim", "delimi", "delimit", "delimite") local subcmd "delimited"
    if inlist("`subcmd'", "exc", "exce") local subcmd "excel"

    if "`subcmd'" == "sasxport" & _caller()<16 local subcmd "sasxport5"
    if inlist("`subcmd'","hav","have") local subcmd "haver"
    if inlist("`subcmd'","haverd","haverdi","haverdir","haverdire","haverdirec") local subcmd "haverdirect"
    if inlist("`subcmd'","haver","haverdirect") & "`c(os)'"!="Windows" {
        di as error "import `subcmd' is not supported on this platform."
        exit 198
    }

    if inlist("`subcmd'","fred","haver","haverdirect") {
        di as error "cimport `subcmd': a C provider adapter is not yet available"
        exit 198
    }
    if inlist("`subcmd'", "sas", "spss", "sasxport5", "sasxport8", "dbase", "shp") {
        _cio_import `subcmd' `0'
        return add
        exit
    }

    if "`subcmd'" == "excel" {
        * Dispatch to Excel import handler
        cimport_excel `0'
        return add
        exit
    }
    else if "`subcmd'" != "delimited" {
        di as error "cimport: unknown subcommand `subcmd'"
        di as error "Supported file formats: delimited, excel, sas, spss, sasxport5, sasxport8, dbase, shp"
        exit 198
    }

    * Now parse the rest with import delimited-style syntax
    * Support both "using filename" and just "filename" (like import delimited does)
    syntax [anything] [using/] [, Delimiters(string asis) VARNames(string) CLEAR ///
        CASE(string) ENCoding(string) BINDQuotes(string) ///
        STRIPQuotes(string) ROWRange(string) COLRange(string) Verbose THReads(integer 0) ///
        ASFloat ASDOUBle NUMERICcols(string) STRINGcols(string) ///
        DECIMALSEParator(string) GROUPSEParator(string) ///
        MAXQUOTEDrows(string) EMPTYlines(string) ///
        PARSELocale(string) CHARSET(string) FAVORSTRFixed COLLAPSEDelimiters ///
        SNAPshot FAST]
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
    if _rc exit 198
    if `maxquotedrows'<0 exit 198
    foreach opt in numericcols stringcols {
        if "``opt''" != "" & "``opt''" != "_all" {
            numlist "``opt''", integer range(>0) sort
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

    * Native delimited import supplies .csv when no extension is given.
    mata: st_local("__basename", pathbasename(st_local("using")))
    if !strpos(`"`__basename'"', ".") local using `"`using'.csv"'
    * Validate file exists
    confirm file `"`using'"'

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
        _cimport_range `"`rowrange'"'
        local startrow = r(first)
        local endrow = r(last)
    }

    local __cimport_startrow `startrow'
    local __cimport_endrow `endrow'
    * Parse colrange option
    local startcol = 0
    local endcol = 0
    if "`colrange'" != "" {
        _cimport_range `"`colrange'"'
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

    * Parse emptylines option - default is skip
    if "`emptylines'" == "include" local emptylines "fill"
    if "`emptylines'" != "" {
        if !inlist("`emptylines'", "skip", "fill") {
            di as error "cimport: emptylines() must be skip or fill"
            exit 198
        }
    }
    else {
        local emptylines "skip"
    }

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

    * Record start time
    local __ct_elapsed99 = 0
    if "`__ct_elapsed99'" == "" local __ct_elapsed99 = 0
    quietly _ctools_clock __ct_start99


    * =========================================================================
    * PHASE 1: Scan the CSV to get metadata
    * =========================================================================

    if "`__ct_elapsed11'" == "" local __ct_elapsed11 = 0

    quietly _ctools_clock __ct_start11

    local __cimport_filename `"`using'"'
    capture noisily plugin call ctools_plugin, ///
        "cimport `threads_code' scan @filename `plugin_delim' `opt_noheader' `opt_headerrow' `opt_verbose' `opt_bindquotes' `opt_asfloat' `opt_asdouble' `opt_decimalsep' `opt_groupsep' `opt_emptylines' `opt_maxquotedrows' `encoding_opt'"

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
    local longtypes ${_cimport_longtypes}
    local hasheader = 0${_cimport_hasheader}
    mata: st_local("returned_encoding",st_global("_cimport_encoding"))
    mata: st_local("returned_delimiters",subinstr(st_global("_cimport_delimiters"),char(9),char(92)+"t"))

    * Clean up global macros
    macro drop _cimport_nobs _cimport_nvar _cimport_varnames ///
               _cimport_vartypes _cimport_numtypes _cimport_strlens _cimport_longtypes ///
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
        exit 103
    }
    if `nobs' > 2147483647 {
        di as error "ctools does not support datasets exceeding 2^31 (2.147 billion) observations"
        di as error "This is a limitation of Stata's plugin API"
        capture plugin call ctools_plugin, "cimport clear"
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
            exit `rc'
        }
    }
    if "`clear'" != "" clear

    * Handle empty file or file with header only - match Stata's behavior (rc=0, N=0, k=0)
    if `nvar' == 0 {
        capture plugin call ctools_plugin, "cimport clear"
        quietly _ctools_clock __ct_now
        local __ct_elapsed11 = `__ct_elapsed11' + (`__ct_now' - `__ct_start11')
        quietly _ctools_clock __ct_now
        local __ct_elapsed99 = `__ct_elapsed99' + (`__ct_now' - `__ct_start99')
        * Elapsed seconds are held in caller-local timing state.
        local elapsed = `__ct_elapsed99'
        di as text "(0 vars, 0 obs)"
        local __ct_elapsed11 = 0
        local __ct_elapsed99 = 0
        return scalar N = 0
        return scalar k = 0
        return scalar time = `elapsed'
        return local filename `"`using'"'
        return local encoding `"`macval(returned_encoding)'"'
        return local delimiters `"`macval(returned_delimiters)'"'
        if "`__cimport_frame'" != "" frame drop `__cimport_frame'
        else if `__cimport_snap' restore, not
        exit 0
    }
    quietly _ctools_clock __ct_now
    local __ct_elapsed11 = `__ct_elapsed11' + (`__ct_now' - `__ct_start11')

    * Phases 2 and 3 run as one block: a failure or Break in them restores
    * the data in memory
    capture noisily {


    * =========================================================================
    * PHASE 2: Create variables in Stata
    * =========================================================================

    if "`__ct_elapsed12'" == "" local __ct_elapsed12 = 0

    quietly _ctools_clock __ct_start12

    * Set number of observations as 1 to create empty variables, change after
    quietly set obs 1

    * Create variables with appropriate types
    * Numeric subtypes: 0=double, 1=float, 2=long, 3=int, 4=byte
    local i = 1
    foreach vname of local varnames {
        * Taken from the front of each list: word `i' of costs O(i) per column.
        gettoken vtype vartypes : vartypes
        gettoken ntype numtypes : numtypes
        gettoken vlen strlens : strlens
        gettoken vlong longtypes : longtypes

        * The plugin names ASCII headers. A header with other characters is
        * named in Mata as native does: Unicode letters, decimal digits (not
        * leading) and underscores, cut to 32 characters.
        mata: st_local("vname", _cimport_name(st_local("vname"), st_global("_cimport_label_`i'")))
        if "`vname'" == "" local vname "v`i'"
        * Apply case transformation (after the cut, as native)
        local vprefix "v"
        if "`case'" == "lower" {
            local vname = ustrlower("`vname'")
        }
        else if "`case'" == "upper" {
            local vname = ustrupper("`vname'")
            local vprefix "V"
        }

        * The data hold only the variables named so far, so a name is taken
        * when a variable has exactly that name (a name list costs O(k) each).
        mata: st_local("taken", strofreal(_st_varindex(st_local("vname")) < .))
        if `taken' {
            local vname "`vprefix'`i'"
            local suffix = 0
            mata: st_local("taken", strofreal(_st_varindex(st_local("vname")) < .))
            while `taken' {
                local ++suffix
                local vname "`vprefix'`i'_`suffix'"
                mata: st_local("taken", strofreal(_st_varindex(st_local("vname")) < .))
            }
        }

        local storage "double"
        if `vtype' == 1 {
            local storage "str`=max(1,`vlen')'"
            if `vlong' {
                local storage "strL"
                local longcolumns `longcolumns' `i'
            }
        }
        else if `ntype' == 4 local storage "byte"
        else if `ntype' == 3 local storage "int"
        else if `ntype' == 2 local storage "long"
        else if `ntype' == 1 local storage "float"

        * _all is legal in imported data but gen treats it as a varlist token.
        if "`vname'" == "_all" {
            quietly mata: st_addvar(st_local("storage"), st_local("vname"))
        }
        else {
            if `vtype' == 1 capture quietly gen `storage' `vname' = ""
            else capture quietly gen `storage' `vname' = .
            if _rc {
                local vname "`vprefix'`i'"
                local suffix = 0
                mata: st_local("taken", strofreal(_st_varindex(st_local("vname")) < .))
                while `taken' {
                    local ++suffix
                    local vname "`vprefix'`i'_`suffix'"
                    mata: st_local("taken", strofreal(_st_varindex(st_local("vname")) < .))
                }
                if `vtype' == 1 quietly gen `storage' `vname' = ""
                else quietly gen `storage' `vname' = .
            }
        }

        * Header text comes from the file: it is compared and stored in Mata
        * only, never pasted into a Stata expression or command line.
        * (A new variable has no label, so "" leaves it as it is.)
        mata: st_varlabel(st_local("vname"), st_global("_cimport_label_`i'") != st_local("vname") ? st_global("_cimport_label_`i'") : "")
        capture macro drop _cimport_label_`i'
        if `vtype' == 1 format `vname' %`=cond(`vlong',9,max(9,`vlen'))'s
        else {
            local storage : type `vname'
            if inlist("`storage'", "byte", "int") format `vname' %8.0g
            else if "`storage'" == "long" format `vname' %12.0g
            else if "`storage'" == "double" format `vname' %10.0g
        }
        local i = `i' + 1
    }

    * Set number of observations
    if `nobs' quietly set obs `nobs'
    else quietly drop in 1

    quietly _ctools_clock __ct_now

    local __ct_elapsed12 = `__ct_elapsed12' + (`__ct_now' - `__ct_start12')

    * =========================================================================
    * PHASE 3: Load data into variables
    * =========================================================================

    if "`__ct_elapsed13'" == "" local __ct_elapsed13 = 0

    quietly _ctools_clock __ct_start13

    * Get list of all variables we just created (in order)
    unab allvars : *

    * Each strL column arrives in one plugin call, through a file of
    * length-prefixed values read in Mata.
    if `nobs' & "`longcolumns'" != "" {
        tempfile __cimport_blobpath
        foreach column of local longcolumns {
            plugin call ctools_plugin, "cimport blobfile `column'"
            mata: _cimport_read_blobs(st_local("__cimport_blobpath"), `column', `nobs')
        }
        capture erase `"`__cimport_blobpath'"'
    }

    if !`nobs' plugin call ctools_plugin, "cimport clear"
    if `nobs' capture noisily plugin call ctools_plugin *, ///
        "cimport `threads_code' load @filename `plugin_delim' `opt_noheader' `opt_headerrow' `opt_verbose' `opt_bindquotes' `opt_asfloat' `opt_asdouble' `opt_decimalsep' `opt_groupsep' `opt_emptylines' `opt_maxquotedrows' `encoding_opt'"

    local load_rc = _rc
    if `load_rc' {
        di as error "Error loading CSV data (rc=`load_rc'); imported data are incomplete"
        exit `load_rc'
    }

    quietly _ctools_clock __ct_now

    local __ct_elapsed13 = `__ct_elapsed13' + (`__ct_now' - `__ct_start13')

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
        nobreak {
            if "`__cimport_frame'" != "" frame copy `__cimport_frame' `c(frame)', replace
            else if `__cimport_snap' restore
        }
        exit `rc'
    }

    quietly _ctools_clock __ct_now

    local __ct_elapsed99 = `__ct_elapsed99' + (`__ct_now' - `__ct_start99')
    * Elapsed seconds are held in caller-local timing state.
    local elapsed = `__ct_elapsed99'

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
        capture scalar drop _cimport_time_mmap _cimport_time_encoding _cimport_time_parse _cimport_time_infer
        capture scalar drop _cimport_time_cache _cimport_time_store _cimport_time_total
        capture scalar drop _cimport_threads_max _cimport_openmp_enabled
    }

    local __ct_elapsed11 = 0
    local __ct_elapsed12 = 0
    local __ct_elapsed13 = 0
    local __ct_elapsed99 = 0

    if "`__cimport_frame'" != "" frame drop `__cimport_frame'
    else if `__cimport_snap' restore, not

    * Return results
    return scalar N = _N
    return scalar k = c(k)
    return scalar time = `elapsed'
    return local filename `"`using'"'
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
    syntax [anything(equalok)] [using/] [, SHEET(string) CELLRange(string) FIRSTrow ///
        ALLSTRING(string) ALLString CASE(string) CLEAR Verbose DEScribe LOCALE(string) DETAIL ///
        SNAPshot FAST]

    * Native describe is an exclusive workbook-inspection syntax.
    if "`describe'"!="" & (`"`sheet'`cellrange'`firstrow'`allstring'`case'`clear'`locale'`detail'"'!="") exit 198

    if `"`allfmt'"' != "" {
        if "`allstring'" == "allstring" exit 198
        mata: st_local("fmtvalid",strofreal(st_isnumfmt(st_local("allfmt"))))
        if !`fmtvalid' exit 198
        local allstring "allstring"
    }

    if substr("`cellrange'",-1,1)==":" exit 198

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
    if "`describe'"!="" & `"`extvarlist'"'!="" exit 198

    * With no extension, native import tries .xls before .xlsx.
    mata: st_local("base",pathbasename(st_local("using")))
    if !strpos(`"`base'"',".") {
        capture confirm file `"`using'.xls"'
        if !_rc local using `"`using'.xls"'
        else local using `"`using'.xlsx"'
    }
    confirm file `"`using'"'
    if `"`extvarlist'"' != "" & ("`case'" != "" | "`firstrow'" != "") exit 198

    if lower(substr(`"`using'"',-4,4)) == ".xls" {
        local opts ""
        if `"`allfmt'"' != "" local opts `"allstringfmt(`"`allfmt'"')"'
        if `"`sheet'"' != "" local opts `"`opts' sheet(`"`sheet'"')"'
        if "`cellrange'" != "" local opts "`opts' cellrange(`cellrange')"
        if "`case'" != "" local opts "`opts' case(`case')"
        _cio_import xls `extvarlist' using `"`using'"', `firstrow' `allstring' `clear' `describe' `snapshot' `fast' `opts'
        return add
        exit
    }

    * Validate .xlsx extension
    local ext = substr(`"`using'"', -5, 5)
    if lower("`ext'") != ".xlsx" {
        di as error "cimport excel: file must have .xlsx extension"
        exit 198
    }

    local opts ""
        if `"`allfmt'"' != "" local opts `"allstringfmt(`"`allfmt'"')"'
    if `"`sheet'"' != "" local opts `"`opts' sheet(`"`sheet'"')"'
    if "`cellrange'" != "" local opts "`opts' cellrange(`cellrange')"
    if "`case'" != "" local opts "`opts' case(`case')"
    _cio_import xlsx `extvarlist' using `"`using'"', `firstrow' `allstring' `clear' `describe' `snapshot' `fast' `opts'
    return add
    exit

end

program define _cimport_delimiters, rclass
    version 14.1
    syntax anything(name=delimiters) [, COLLApse ASSTRing]
    if "`collapse'"!="" & "`asstring'"!="" exit 198
    local quoted = strpos(`"`macval(delimiters)'"',char(34))>0
    gettoken delimiters extra : delimiters
    if strtrim(`"`extra'"')!="" | `"`delimiters'"'=="" exit 198
    local delimiters = subinstr(`"`macval(delimiters)'"',"\t",char(9),.)
    if !`quoted' {
        if !inlist("`delimiters'","tab","comma","space","whitespace",char(9),",") exit 198
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
/* cname is the plugin's name for an ASCII header. Other headers keep only
   Unicode letters, decimal digits and underscores, drop leading digits and
   are cut to 32 characters, as native import delimited names them. */
string scalar _cimport_name(string scalar cname, string scalar header) {
    string scalar name
    if (ustrinvalidcnt(cname) == 0 & ustrlen(cname) == strlen(cname)) return(cname)
    name = ustrregexrf(ustrregexra(header, "[^\p{L}\p{Nd}_]", ""), "^\p{Nd}+", "")
    return(usubstr(name, 1, 32))
}
/* Native infers display formats from a sample of the file's lines: the
   first 513 (header included), then the later lines in groups of 512 counted
   from the second line, group b when b is a power of two from 2, and a last
   group that is not full. Integers there display as %8.0g, or %12.0g in long
   range; a non-integer as float (%9.0g), or double (%10.0g) with asdouble or
   an integer beyond long range in the sample. */
void _cimport_sample_formats(real scalar hasheader, real scalar asdouble) {
    real scalar j, n, b, last, lo, hi
    real colvector obs, x, y
    n = st_nobs()
    if (n < 1) return
    obs = 1::min((n, 513 - hasheader))
    for (b = 2; 512 * b + 2 - hasheader <= n; b = 2 * b) {
        obs = obs \ ((512 * b + 2 - hasheader)::min((n, 512 * b + 513 - hasheader)))
    }
    last = n - 2 + hasheader
    if (last >= 512 & mod(last + 1, 512)) obs = obs \ ((512 * floor(last / 512) + 2 - hasheader)::n)
    obs = uniqrows(obs)
    for (j = 1; j <= st_nvar(); j++) {
        if (!st_isnumvar(j)) continue
        x = st_data(obs, j)
        x = select(x, x :< .)
        if (!length(x)) continue
        y = select(x, x :== floor(x))
        lo = length(y) ? min(y) : 0
        hi = length(y) ? max(y) : 0
        if (any(x :!= floor(x))) st_varformat(j, asdouble | lo < -2147483647 | hi > 2147483620 ? "%10.0g" : "%9.0g")
        else if (lo >= -32767 & hi <= 32740) st_varformat(j, "%8.0g")
        else if (lo >= -2147483647 & hi <= 2147483620) st_varformat(j, "%12.0g")
        else st_varformat(j, "%10.0g")
    }
}
void _cimport_read_blobs(string scalar path, real scalar col, real scalar n) {
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

* Native row/column ranges support f:l as well as omitted bounds.
program define _cimport_range, rclass
    version 14.1
    args range
    local colon = strpos("`range'",":")
    local first 0
    local last 0
    if `colon' {
        local left = strtrim(substr("`range'",1,`colon'-1))
        local right = strtrim(substr("`range'",`colon'+1,.))
        if "`left'"=="f" local left "1"
        if "`left'"!="" {
            local first = real("`left'")
            if missing(`first') | `first'<1 | strpos("`left'",".") exit 198
        }
        if inlist("`right'","","l") & "`left'"!="" local last 0
        else {
            local last = real("`right'")
            if missing(`last') | `last'<1 | (`colon'==1 & strpos("`right'",".")) exit 198
        }
        if `last' & `first'>`last' exit 198
    }
    else {
        local first = real("`range'")
        if missing(`first') | `first'<1 | strpos("`range'",".") exit 198
    }
    return scalar first=`first'
    return scalar last=`last'
end
