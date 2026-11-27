*! version 1.4.0 20261127 github.com/mdroste/stata-ctools
*! C-accelerated fuzzy joins: record linkage, similarity joins and grouping of similar strings (ctools)
*! original: reclink 1.7 (January 14, 2010) by Michael Blasnik
*! original: matchit 1.5.2 (May 2020) by Julio D. Raffo
*! original: strgroup 1.0.4 (December 17, 2021) by Julian Reif
*! original: see help cfuzzyjoin for details

* cfuzzyjoin link   reclink's record linkage: the best fuzzy match of each
*                   master observation, scored from several variables
* cfuzzyjoin match  matchit's similarity join of two datasets (or the
*                   similarity of two columns): every pair above a threshold
* cfuzzyjoin group  strgroup's groups of similar strings (Levenshtein)
*
* The ado file parses, checks and lays out the data; finding the matches,
* scoring them and writing the results is done by the plugin.

capture mata: mata drop _cfj_strl_hold()
capture mata: mata drop _cfj_strl_fill()
mata:
/* Hold the strL variables of the using data (all observations) */
void _cfj_strl_hold(string scalar vars)
{
    external string matrix _cfj_strl_values
    _cfj_strl_values = st_sdata(., tokens(vars))
}

/* Copy them into the output: row r takes using observation urow[r] */
void _cfj_strl_fill(string scalar vars, string scalar urow)
{
    external string matrix _cfj_strl_values
    real colvector u
    real scalar i, j
    string rowvector v

    v = tokens(vars)
    u = st_data(., urow)
    for (i = 1; i <= rows(u); i++) {
        if (u[i] >= .) continue
        for (j = 1; j <= cols(v); j++) st_sstore(i, v[j], _cfj_strl_values[u[i], j])
    }
    _cfj_strl_values = J(0, 0, "")
}
end

program define cfuzzyjoin, rclass byable(onecall)
    version 14.1
    if _N > 2147483647 {
        di as error "ctools does not support datasets exceeding 2^31 (2.147 billion) observations"
        di as error "This is a limitation of Stata's plugin API"
        exit 920
    }
    gettoken sub rest : 0, parse(" ,")
    local byvars
    if _by() local byvars `_byvars'
    if "`sub'" == "group" {
        _cfuzzyjoin_group `rest' byvars(`byvars')
    }
    else if inlist("`sub'", "link", "match") {
        if "`byvars'" != "" {
            di as error "cfuzzyjoin `sub' may not be combined with by"
            exit 190
        }
        _cfuzzyjoin_`sub' `rest'
    }
    else {
        di as error `"cfuzzyjoin: unknown subcommand `sub'; use link, match or group"'
        exit 198
    }
    return add
end

* -----------------------------------------------------------------------------
* Loading the plugin
* -----------------------------------------------------------------------------
program define _cfuzzyjoin_plugin
    _ctools_load
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    c_local __cfj_plugin "`__ctools_plugin'"
end

* =============================================================================
* group: strgroup
* =============================================================================
program define _cfuzzyjoin_group, rclass
    version 14.1
    syntax varname [if] [in], GENerate(name) THRESHold(real) [first NORMalize(string) ///
        noCLEAN force Verbose THReads(integer 0) byvars(string)]
    if !inrange(`threads', 0, 256) {
        di as error "threads() must be between 0 and 256"
        exit 198
    }
    quietly _ctools_clock __ct_start
    * strgroup excludes missing strings (they would spoil the normalization)
    marksample touse, strok
    confirm new var `generate'
    confirm string var `varlist'

    local normalize `"`=lower(`"`normalize'"')'"'
    if !inlist(`"`normalize'"', "", "longer", "shorter", "none") {
        di as error `"Normalization '`normalize'' not allowed"'
        exit 198
    }
    * strgroup means to require an integer threshold under normalize(none),
    * but its check runs after normalize() has become a number and never
    * fires: a fractional threshold is accepted (a distance of at most 1.5 is
    * a distance of at most 1)
    if `threshold' <= 0 {
        di as error "{cmd:threshold(}{it:#}{cmd:)} must be greater than 0"
        exit 198
    }
    local norm = cond("`normalize'" == "longer", 1, cond("`normalize'" == "none", 2, 0))
    * by-groups must be in sort order, as by requires
    local nby : word count `byvars'
    if `nby' {
        local sorted : sortedby
        local lead
        forvalues j = 1/`nby' {
            local lead `lead' `: word `j' of `sorted''
        }
        if "`: list retokenize lead'" != "`: list retokenize byvars'" {
            di as error "not sorted"
            exit 5
        }
    }
    local threadopt
    if `threads' > 0 local threadopt threads(`threads')
    local flags
    if "`first'" != "" local flags `flags' first
    if "`clean'" != "noclean" local flags `flags' clean
    if "`verbose'" != "" local flags `flags' verbose

    _cfuzzyjoin_plugin
    quietly _ctools_lap __time_setup __ct_start `__ct_start'
    quietly gen long `generate' = .
    _ctools_strw `varlist' `byvars' `generate' `touse'
    capture noisily plugin call ctools_plugin `varlist' `byvars' `generate' `touse', ///
        "cfuzzyjoin `threadopt' group nby=`nby' normalize=`norm' threshold(`threshold') `flags'"
    local rc = _rc
    if `rc' {
        capture drop `generate'
        exit `rc'
    }
    quietly _ctools_lap __time_plugin __ct_start `__ct_start'
    return scalar N = scalar(_cfuzzyjoin_N)
    return scalar groups = scalar(_cfuzzyjoin_ngroups)
    return scalar distinct = scalar(_cfuzzyjoin_nunique)
    if "`verbose'" != "" {
        di as text ""
        di as text "{hline 60}"
        di as text "cfuzzyjoin group timing breakdown:"
        di as text "{hline 60}"
        di as text "  Stata setup:              " as result %8.4f `__time_setup' " sec"
        di as text "  C plugin phases:"
        di as text "    Data load:              " as result %8.4f scalar(_cfuzzyjoin_t_load) " sec"
        di as text "    Grouping:               " as result %8.4f scalar(_cfuzzyjoin_t_compute) " sec"
        di as text "    Data store:             " as result %8.4f scalar(_cfuzzyjoin_t_store) " sec"
        di as text "  Plugin call:              " as result %8.4f `__time_plugin' " sec"
        di as text "{hline 60}"
        di as text "  Strings:          " as result %12.0fc scalar(_cfuzzyjoin_N)
        di as text "  Distinct strings: " as result %12.0fc scalar(_cfuzzyjoin_nunique)
        di as text "  Distances:        " as result %12.0fc scalar(_cfuzzyjoin_pairs)
        di as text "  Groups:           " as result %12.0fc scalar(_cfuzzyjoin_ngroups)
        di as text "  Threads:          " as result %12.0f scalar(_cfuzzyjoin_threads)
        di as text "{hline 60}"
    }
    capture _ctools_cleanup cfuzzyjoin
end

* =============================================================================
* match: matchit
* =============================================================================
program define _cfuzzyjoin_match, rclass
    version 14.1
    syntax varlist(min=2 max=2) [using/] [, IDUsing(name) TXTUsing(name) ///
        SIMilmethod(string asis) Weights(string) WGTFile(string) Score(string) ///
        Threshold(real .5) Flag(real 20) DIagnose STOPWordsauto SWThreshold(real .2) ///
        OVERride Generate(string) KEEPMata TIme Verbose THReads(integer 0)]
    if !inrange(`threads', 0, 256) {
        di as error "threads() must be between 0 and 256"
        exit 198
    }
    quietly _ctools_clock __ct_start
    if "`using'" == "" {
        local match columns
        tokenize `varlist'
        local str1 `1'
        local str2 `2'
        di "Matching columns `str1' and `str2'"
        confirm string variable `str1'
        confirm string variable `str2'
    }
    else {
        local match index
        tokenize `varlist'
        local idmaster `1'
        local txtmaster `2'
        di "Matching current dataset with `using'"
        confirm numeric variable `idmaster'
        confirm string variable `txtmaster'
        local curtype1 : type `txtmaster'
        * use's rule: .dta is assumed when no extension is given
        local ufile `"`using'"'
        capture confirm file `"`ufile'"'
        if _rc & !strpos(substr(`"`ufile'"', length(`"`ufile'"') - 3, 4), ".") local ufile `"`ufile'.dta"'
        confirm file `"`ufile'"'
        if "`idusing'" == "" | "`txtusing'" == "" {
            di as error "options idusing() and txtusing() are required with using"
            exit 198
        }
        quietly describe
        if r(changed) > 0 & "`override'" == "" {
            di " "
            di "(!) Unsaved changes will be destroyed after running matching procedure."
            di "    (note: use OVERRIDE option to bypass warning)"
            error 4
        }
    }
    if "`weights'" == "" local weights noweights
    local wopt none
    if "`weights'" != "noweights" {
        if !inlist("`weights'", "simple", "root", "log") {
            di "`weights' not found as a weights function. Check spelling."
            exit 3499
        }
        di "Applying weights function: `weights'"
        local wopt `weights'
        if `"`wgtfile'"' != "" {
            local wfile `"`wgtfile'"'
            capture confirm file `"`wfile'"'
            if _rc & !strpos(substr(`"`wfile'"', length(`"`wfile'"') - 3, 4), ".") local wfile `"`wfile'.dta"'
            confirm file `"`wfile'"'
        }
    }
    else local wgtfile
    * similmethod(): a name, then the function's arguments after a comma
    local __cfj_sim `"`similmethod'"'
    if `"`__cfj_sim'"' == "" local __cfj_sim bigram
    gettoken similfunc : __cfj_sim, parse(", ")
    di "Similarity function: `similfunc'"
    if "`score'" == "" local score jaccard
    if !inlist("`score'", "jaccard", "simple", "minsimple") {
        di "`score' not found as a score computing function. Check spelling."
        exit 3499
    }
    if "`time'" != "" di c(current_date) " " c(current_time)
    if `flag' <= 0 {
        di "Flag cannot equal or lower than zero."
        exit
    }
    local threadopt
    if `threads' > 0 local threadopt threads(`threads')
    local common score(`score') weights(`wopt')
    if `"`wgtfile'"' != "" local common `common' wgtfile
    if "`verbose'" != "" local common `common' verbose

    _cfuzzyjoin_plugin
    quietly _ctools_lap __time_setup __ct_start `__ct_start'

    * -------------------------------------------------------------------------
    * Two columns
    * -------------------------------------------------------------------------
    if "`match'" == "columns" {
        if "`generate'" == "" local similscore similscore
        else local similscore `generate'
        capture confirm new variable `similscore'
        if _rc {
            di "`similscore' does not seem a valid name for the score variable."
            error _rc
        }
        if `"`wgtfile'"' != "" {
            di "Loading `wgtfile' as weights file..." _continue
            preserve
            quietly use `"`wfile'"', clear
            confirm string variable grams
            confirm numeric variable freq
            _ctools_strw grams freq
            plugin call ctools_plugin grams freq, "cfuzzyjoin `threadopt' match_weights"
            restore
            di "Done!"
        }
        quietly gen double `similscore' = .
        _ctools_strw `str1' `str2' `similscore'
        capture noisily plugin call ctools_plugin `str1' `str2' `similscore', ///
            "cfuzzyjoin `threadopt' match_columns `common'"
        local rc = _rc
        if `rc' {
            capture drop `similscore'
            capture plugin call ctools_plugin, "cfuzzyjoin clear"
            exit `rc'
        }
        di "Done!"
        quietly _ctools_lap __time_plugin __ct_start `__ct_start'
        if "`time'" != "" di c(current_date) " " c(current_time)
        if "`verbose'" != "" {
            di as text ""
            di as text "{hline 60}"
            di as text "cfuzzyjoin match timing breakdown:"
            di as text "{hline 60}"
            di as text "  Stata setup:              " as result %8.4f `__time_setup' " sec"
            di as text "    Grams:                  " as result %8.4f scalar(_cfuzzyjoin_t_grams) " sec"
            di as text "    Scores and store:       " as result %8.4f scalar(_cfuzzyjoin_t_compute) " sec"
            di as text "  Plugin calls:             " as result %8.4f `__time_plugin' " sec"
            di as text "  Threads:      " as result %12.0f scalar(_cfuzzyjoin_threads)
            di as text "{hline 60}"
        }
        capture _ctools_cleanup cfuzzyjoin
        exit
    }

    * -------------------------------------------------------------------------
    * Two datasets: the result replaces the data in memory
    * -------------------------------------------------------------------------
    local mopts `common' threshold(`threshold')
    if "`stopwordsauto'" != "" local mopts `mopts' stopwords swthreshold(`swthreshold')
    if "`diagnose'" != "" local mopts `mopts' diagnose
    preserve
    _ctools_strw `idmaster' `txtmaster'
    plugin call ctools_plugin `idmaster' `txtmaster', "cfuzzyjoin `threadopt' match_master"
    if `"`wgtfile'"' != "" {
        di "Loading `wgtfile' as weights file..." _continue
        quietly use `"`wfile'"', clear
        confirm string variable grams
        confirm numeric variable freq
        _ctools_strw grams freq
        plugin call ctools_plugin grams freq, "cfuzzyjoin `threadopt' match_mweights"
        di "Done!"
    }
    di "Loading USING file: `using'"
    quietly use `"`ufile'"', clear
    confirm numeric variable `idusing'
    confirm string variable `txtusing'
    local curtype2 : type `txtusing'
    if "`stopwordsauto'" != "" {
        di "Generating stopwords automatically, threshold set at:" `swthreshold'
    }
    di "Indexing USING file."
    di "Computing results"
    quietly _ctools_lap __time_load __ct_start `__ct_start'
    _ctools_strw `idusing' `txtusing'
    capture noisily plugin call ctools_plugin `idusing' `txtusing', ///
        "cfuzzyjoin `threadopt' match_using `mopts'"
    local rc = _rc
    if `rc' {
        capture plugin call ctools_plugin, "cfuzzyjoin clear"
        exit `rc'
    }
    if "`stopwordsauto'" != "" di "Done!"
    quietly _ctools_lap __time_compute __ct_start `__ct_start'

    * Names: matchit renames the using variables (and the score) that clash
    * with the master's
    if "`generate'" == "" local similscore similscore
    else local similscore `generate'
    local i 1
    local vartemp `idusing'
    while inlist("`vartemp'", "`idmaster'", "`txtmaster'") {
        local vartemp `idusing'`i'
        local ++i
    }
    local idusing `vartemp'
    local i 1
    local vartemp `txtusing'
    while inlist("`vartemp'", "`idmaster'", "`txtmaster'", "`idusing'") {
        local vartemp `txtusing'`i'
        local ++i
    }
    local txtusing `vartemp'
    local i 1
    local vartemp `similscore'
    while inlist("`vartemp'", "`idmaster'", "`txtmaster'", "`idusing'", "`txtusing'") {
        local vartemp `similscore'`i'
        local ++i
    }
    local similscore `vartemp'

    * Storage types: those compress leaves (matchit compresses its result)
    local npairs = scalar(_cfuzzyjoin_npairs)
    local tnames byte int long float double
    if `npairs' > 0 {
        local t1 : word `=scalar(_cfuzzyjoin_type_idm)' of `tnames'
        local t2 = "str" + string(scalar(_cfuzzyjoin_len_txtm))
        local t3 : word `=scalar(_cfuzzyjoin_type_idu)' of `tnames'
        local t4 = "str" + string(scalar(_cfuzzyjoin_len_txtu))
        local t5 : word `=scalar(_cfuzzyjoin_type_score)' of `tnames'
        if scalar(_cfuzzyjoin_len_txtm) > 2045 local t2 strL
        if scalar(_cfuzzyjoin_len_txtu) > 2045 local t4 strL
    }
    else {
        local t1 double
        local t2 `curtype1'
        local t3 double
        local t4 `curtype2'
        local t5 double
    }
    nobreak {
        quietly drop _all
        label drop _all
        quietly set obs `npairs'
        mata: (void) st_addvar(tokens("`t1' `t2' `t3' `t4' `t5'"), ///
            tokens("`idmaster' `txtmaster' `idusing' `txtusing' `similscore'"), 1)
        * matchit makes them double, and compress keeps the format
        format `idmaster' `idusing' `similscore' %10.0g
        capture noisily plugin call ctools_plugin `idmaster' `txtmaster' `idusing' `txtusing' `similscore', ///
            "cfuzzyjoin `threadopt' match_store"
        local rc = _rc
    }
    if `rc' {
        capture plugin call ctools_plugin, "cfuzzyjoin clear"
        exit `rc'
    }
    di "{tab}Done!"
    di "Total search space saved by index: " scalar(_cfuzzyjoin_saved) "%"
    restore, not
    global S_FN
    global S_FNDATE
    quietly _ctools_lap __time_store __ct_start `__ct_start'
    if "`time'" != "" di c(current_date) " " c(current_time)
    return scalar N = `npairs'
    if "`verbose'" != "" {
        di as text ""
        di as text "{hline 60}"
        di as text "cfuzzyjoin match timing breakdown:"
        di as text "{hline 60}"
        di as text "  Stata setup:              " as result %8.4f `__time_setup' " sec"
        di as text "  Load master and using:    " as result %8.4f `__time_load' " sec"
        di as text "  C plugin phases:"
        di as text "    Grams:                  " as result %8.4f scalar(_cfuzzyjoin_t_grams) " sec"
        di as text "    Index:                  " as result %8.4f scalar(_cfuzzyjoin_t_index) " sec"
        di as text "    Scores:                 " as result %8.4f scalar(_cfuzzyjoin_t_compute) " sec"
        di as text "  Score call:               " as result %8.4f `__time_compute' " sec"
        di as text "  Result dataset:           " as result %8.4f `__time_store' " sec"
        di as text "{hline 60}"
        di as text "  Distinct grams:   " as result %12.0fc scalar(_cfuzzyjoin_ngrams)
        di as text "  Stopwords:        " as result %12.0fc scalar(_cfuzzyjoin_nstop)
        di as text "  Pairs scored:     " as result %12.0fc scalar(_cfuzzyjoin_pairs)
        di as text "  Pairs kept:       " as result %12.0fc `npairs'
        di as text "  Threads:          " as result %12.0f scalar(_cfuzzyjoin_threads)
        di as text "{hline 60}"
    }
    capture _ctools_cleanup cfuzzyjoin
end

* =============================================================================
* link: reclink
* =============================================================================
program define _cfuzzyjoin_link, rclass
    version 14.1
    syntax varlist using/ , IDMaster(str) IDUsing(str) Gen(str) [WMatch(str) WNOMatch(str) ///
        EXClude(str) UVarlist(str) MINScore(real 0.6) MINBigram(real 0.6) ORBLock(str) ///
        REQuired(str) _Merge(str) UPRefix(str) EXActstr(str) DEBug(str) ///
        Verbose THReads(integer 0)]
    if !inrange(`threads', 0, 256) {
        di as error "threads() must be between 0 and 256"
        exit 198
    }
    quietly _ctools_clock __ct_start
    * reclink runs under version 8.2, whose syntax keeps a real option as text
    * of about seven significant digits: it compares scores with minscore()
    * and scales bigram scores with minbigram() as that text gives them
    _cfuzzyjoin_v8real, minscore(`minscore') minbigram(`minbigram')

    * reclink's checks, in its order
    if "`_merge'" == "" local _merge "_merge"
    capture confirm variable `_merge'
    if _rc == 0 {
        di as err "`_merge' already exists"
        exit 198
    }
    local nvars : word count `varlist'
    if "`wmatch'" == "" local wmatch : di _dup(`nvars') " 1 "
    if "`wnomatch'" == "" local wnomatch "`wmatch'"
    if "`uprefix'" == "" local uprefix "U"
    * (reclink checks the wmatch() weights only)
    foreach val of local wmatch {
        capture confirm number `val'
        if _rc | `val' < 1 {
            di as error " weights must be >=1"
            exit 198
        }
    }
    foreach val of local wnomatch {
        capture confirm number `val'
        if _rc {
            di as error " weights must be numbers"
            exit 198
        }
    }
    local nwms : word count `wmatch'
    local nwnms : word count `wnomatch'
    if `nvars' != `nwnms' | `nvars' != `nwms' {
        di as error "wmatch and wnomatch must each have the same # of elements as varlist"
        exit 198
    }
    if "`uvarlist'" == "" local uvarlist "`varlist'"
    local nuvars : word count `uvarlist'
    if `nvars' != `nuvars' {
        di as error "uvarlist must have the same # of elements as varlist"
        exit 198
    }
    if `: list required in varlist' == 0 {
        di as error "all variables in required() option must also be in the main varlist"
        exit 198
    }
    if "`orblock'" != "none" & `: list orblock in varlist' == 0 {
        di as error "all variables in orblock() option must also be in the main varlist"
        exit 198
    }
    if "`orblock'" == "" & `nvars' > 3 local orblock "`varlist'"
    if "`orblock'" == "none" local orblock
    local orblock : list orblock - required
    confirm variable `idmaster'
    capture confirm variable `gen'
    if !_rc {
        di as error "variable `gen' already defined"
        exit 110
    }
    capture confirm variable `idusing', exact
    if !_rc {
        di as error "variable `idusing' already defined in the master data"
        exit 110
    }
    if "`idmaster'" == "`idusing'" {
        di as error "idmaster() and idusing() must name different variables"
        exit 198
    }

    * Per variable: kind (b bigram, e exact string, n numeric), required,
    * or-blocking, weights
    local __cfj_kind
    local __cfj_req
    local __cfj_orb
    local __cfj_wm
    local __cfj_wnm
    local i 0
    foreach var of local varlist {
        local ++i
        capture confirm string var `var'
        if _rc local kind n
        else if `: list var in exactstr' local kind e
        else local kind b
        local __cfj_kind `__cfj_kind' `kind'
        local __cfj_req `__cfj_req' `: list var in required'
        local __cfj_orb `__cfj_orb' `: list var in orblock'
        local __cfj_wm `__cfj_wm' `: word `i' of `wmatch''
        local __cfj_wnm `__cfj_wnm' `: word `i' of `wnomatch''
    }
    * the scaling of bigram scores, as reclink keeps it in a macro
    local biscale = .5/(1-`minbigram')^(1/3)

    local ufile `"`using'"'
    capture confirm file `"`ufile'"'
    if _rc & !strpos(substr(`"`ufile'"', length(`"`ufile'"') - 3, 4), ".") local ufile `"`ufile'.dta"'
    confirm file `"`ufile'"'
    if `"`exclude'"' != "" {
        local xfile `"`exclude'"'
        capture confirm file `"`xfile'"'
        if _rc & !strpos(substr(`"`xfile'"', length(`"`xfile'"') - 3, 4), ".") local xfile `"`xfile'.dta"'
        confirm file `"`xfile'"'
    }
    local threadopt
    if `threads' > 0 local threadopt threads(`threads')
    local vflag
    if "`verbose'" != "" local vflag verbose

    _cfuzzyjoin_plugin
    quietly _ctools_lap __time_setup __ct_start `__ct_start'
    unab mastervars : _all
    local mobs = _N

    preserve
    _ctools_strw `idmaster' `varlist'
    plugin call ctools_plugin `idmaster' `varlist', "cfuzzyjoin `threadopt' link_master id(`idmaster')"
    if `"`exclude'"' != "" {
        quietly use `idmaster' `idusing' using `"`xfile'"', clear
        _ctools_strw `idmaster' `idusing'
        plugin call ctools_plugin `idmaster' `idusing', "cfuzzyjoin `threadopt' link_exclude"
    }
    quietly use `"`ufile'"', clear
    confirm variable `idusing'
    confirm variable `uvarlist'
    unab usingvars : _all
    * The using variables brought into the master: every one but idusing; a
    * matching variable named as in the master takes the prefix (reclink's
    * rename, which stops with r(110) when the prefixed name is taken); names
    * the master already has keep the master's variable, widened to hold the
    * using type as merge widens it (values unchanged).
    local payload
    local finals
    local strls
    local strlfinal
    local widen
    foreach v of local usingvars {
        if "`v'" == "`idusing'" continue
        local final `v'
        if `: list v in varlist' {
            local pref `uprefix'`v'
            if `: list pref in usingvars' {
                di as error "variable `pref' already defined"
                exit 110
            }
            local final `pref'
        }
        if `: list final in mastervars' | "`final'" == "`gen'" {
            local widen `widen' `final' `: type `v''
            continue
        }
        if "`final'" == "`_merge'" {
            di as error "variable `_merge' already defined in the using data"
            exit 110
        }
        if "`: type `v''" == "strL" {
            local strls `strls' `v'
            local strlfinal `strlfinal' `final'
        }
        else {
            local payload `payload' `v'
            local finals `finals' `final'
        }
    }
    local npay : word count `payload'
    quietly _ctools_lap __time_load __ct_start `__ct_start'
    _ctools_strw `idusing' `uvarlist' `payload'
    capture noisily plugin call ctools_plugin `idusing' `uvarlist' `payload', ///
        "cfuzzyjoin `threadopt' link_using k=`nvars' p=`npay' id(`idusing') minscore(`minscore') minbigram(`minbigram') biscale(`biscale') `vflag'"
    local rc = _rc
    if `rc' {
        capture plugin call ctools_plugin, "cfuzzyjoin clear"
        exit `rc'
    }
    quietly _ctools_lap __time_compute __ct_start `__ct_start'
    local perfectn = scalar(_cfuzzyjoin_perfect)
    local tomatch = scalar(_cfuzzyjoin_tomatch)
    local uobs = scalar(_cfuzzyjoin_using_N)
    local nout = scalar(_cfuzzyjoin_nout)
    local expand = scalar(_cfuzzyjoin_expand)

    * The new variables, with everything that describes them: an empty
    * dataset of idusing and the payload, appended to the master
    tempfile template
    * every value label of the using data goes to the master when the master
    * has none of that name, also those of variables not brought in, as
    * reclink's merge copies them
    quietly label dir
    local ulabels `r(names)'
    local nul 0
    foreach l of local ulabels {
        local ++nul
        tempfile ulab`nul'
        quietly label save `l' using `"`ulab`nul''"'
    }
    local tvars `idusing' `payload' `strls'
    local tfinal `idusing' `finals' `strlfinal'
    quietly keep `tvars'
    order `tvars'
    local j 0
    foreach v of local tvars {
        local ++j
        local f : word `j' of `tfinal'
        if "`f'" != "`v'" rename `v' `f'
    }
    if "`strls'" != "" mata: _cfj_strl_hold("`strlfinal'")
    if _N quietly drop in 1/l
    quietly save `"`template'"', emptyok
    restore, preserve

    tempvar cnt mrow ord urow
    nobreak {
        quietly gen float `gen' = .
        format `gen' %5.4f
        quietly append using `"`template'"'
        while "`widen'" != "" {
            gettoken wv widen : widen
            gettoken wt widen : widen
            _cfuzzyjoin_widen `wv' `wt'
        }
        quietly label dir
        local mlabels `r(names)'
        local j 0
        foreach l of local ulabels {
            local ++j
            if !`: list l in mlabels' quietly do `"`ulab`j''"'
        }
        quietly gen byte `_merge' = .
        quietly gen long `cnt' = .
        plugin call ctools_plugin `cnt', "cfuzzyjoin `threadopt' link_layout"
        local rowvar `cnt'
        if `expand' {
            quietly gen long `mrow' = _n
            quietly expand `cnt'
            local rowvar `mrow'
        }
        quietly gen double `ord' = .
        local uflag
        local uvar
        if "`strls'" != "" {
            quietly gen double `urow' = .
            local uflag urow
            local uvar `urow'
        }
        capture noisily plugin call ctools_plugin `rowvar' `gen' `idusing' `finals' `_merge' `ord' `uvar', ///
            "cfuzzyjoin `threadopt' link_store `uflag'"
        local rc = _rc
        if !`rc' {
            if "`strls'" != "" mata: _cfj_strl_fill("`strlfinal'", "`urow'")
            sort `ord'
            drop `ord' `cnt'
            capture drop `mrow'
            capture drop `urow'
            * each prefixed variable after its master counterpart
            foreach v of local varlist {
                capture confirm variable `uprefix'`v', exact
                if !_rc & `: list v in mastervars' order `uprefix'`v', after(`v')
            }
        }
    }
    if `rc' {
        capture plugin call ctools_plugin, "cfuzzyjoin clear"
        exit `rc'
    }
    restore, not
    quietly _ctools_lap __time_store __ct_start `__ct_start'

    if `"`debug'"' != "" {
        preserve
        quietly keep if `_merge' == 3
        keep `idmaster' `gen' `idusing'
        sort `idmaster'
        quietly save `"`debug'"', replace
        restore
    }
    quietly count if missing(`idusing')
    local nomatch = r(N)
    noi di as result _n "`perfectn' perfect matches found" _n
    if `tomatch' > 1000 {
        noi di as text "Going through `tomatch' observation to assess fuzzy matches, each .=5% complete"
        noi di "...................."
    }
    noi di _n "Added: `idusing'= identifier from `using'   `gen' = matching score"
    noi di as res "Observations:  Master N = `mobs'    `using' N= `uobs' "
    if `"`exclude'"' != "" {
        noi di as res "  # Obs after excluding matches in `exclude':  Master = " ///
            scalar(_cfuzzyjoin_master_left) "    `using': " scalar(_cfuzzyjoin_using_left) " "
    }
    noi di as res "  Unique Master Cases: matched = " `mobs'-`nomatch' " (exact = `perfectn'), unmatched = `nomatch'"
    return scalar N = _N
    return scalar matched = `mobs' - `nomatch'
    return scalar perfect = `perfectn'
    return scalar unmatched = `nomatch'
    if "`verbose'" != "" {
        di as text ""
        di as text "{hline 60}"
        di as text "cfuzzyjoin link timing breakdown:"
        di as text "{hline 60}"
        di as text "  Stata setup:              " as result %8.4f `__time_setup' " sec"
        di as text "  Load master and using:    " as result %8.4f `__time_load' " sec"
        di as text "  C plugin phases:"
        di as text "    Index:                  " as result %8.4f scalar(_cfuzzyjoin_t_index) " sec"
        di as text "    Scores:                 " as result %8.4f scalar(_cfuzzyjoin_t_compute) " sec"
        di as text "    Layout:                 " as result %8.4f scalar(_cfuzzyjoin_t_layout) " sec"
        di as text "  Match call:               " as result %8.4f `__time_compute' " sec"
        di as text "  Output:                   " as result %8.4f `__time_store' " sec"
        di as text "{hline 60}"
        di as text "  Using combinations: " as result %12.0fc scalar(_cfuzzyjoin_ncomb)
        di as text "  Pairs scored:       " as result %12.0fc scalar(_cfuzzyjoin_pairs)
        di as text "  Threads:            " as result %12.0f scalar(_cfuzzyjoin_threads)
        di as text "{hline 60}"
    }
    capture _ctools_cleanup cfuzzyjoin
end

* minscore() and minbigram() as reclink's version 8.2 syntax gives them
program define _cfuzzyjoin_v8real
    version 8.2
    syntax , MINScore(real) MINBigram(real)
    c_local minscore `minscore'
    c_local minbigram `minbigram'
end

* _cfuzzyjoin_widen var type: widen master variable var to hold the values of
* a using variable of the given type, as merge does: the longer str#, strL,
* or the numeric type that holds both (a string and a numeric type, or a
* narrower using type, leave var as it is); values and format are unchanged
program define _cfuzzyjoin_widen
    version 14.1
    args var ut
    local mt : type `var'
    if "`mt'" == "`ut'" exit
    local ms = substr("`mt'", 1, 3) == "str"
    local us = substr("`ut'", 1, 3) == "str"
    if `ms' != `us' | "`mt'" == "strL" exit
    if `ms' {
        if "`ut'" != "strL" & real(substr("`ut'", 4, .)) <= real(substr("`mt'", 4, .)) exit
        local nt `ut'
    }
    else {
        local order byte int long float double
        local a : list posof "`mt'" in order
        local b : list posof "`ut'" in order
        if `a' == 5 | `b' == 5 local nt double
        else if `a' == 4 | `b' == 4 local nt = cond(`a' == 3 | `b' == 3, "double", "float")
        else local nt : word `=max(`a', `b')' of `order'
        if "`nt'" == "`mt'" exit
    }
    quietly recast `nt' `var'
end

* ctools-helpers:begin (generated by scripts/sync_release.py from src/ado; edit those files)
* ctools: drop the scalars and matrices a command and its plugin leave behind
* _ctools_cleanup <command> [<command> ...] drops every scalar and matrix
* named __<command>_* or _<command>_*.
program define _ctools_cleanup
    version 14.1
    foreach command of local 0 {
        foreach pattern in __`command'_* _`command'_* {
            local found : all scalars "`pattern'"
            if "`found'" != "" capture scalar drop `found'
            local found : all matrices "`pattern'"
            if "`found'" != "" capture matrix drop `found'
        }
    }
end

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

* ctools: seconds since the previous lap of a caller-local wall clock
* _ctools_lap <target> <clock> <clock value> stores the elapsed seconds in the
* caller's local <target> and restarts the clock, a local the caller starts
* with _ctools_clock <clock> and passes back expanded:
*     _ctools_clock __ct_lap
*     ...
*     _ctools_lap __time_setup __ct_lap `__ct_lap'
program define _ctools_lap
    version 14.1
    args target clock start
    quietly _ctools_clock now
    c_local `target' = `now' - `start'
    c_local `clock' `now'
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
