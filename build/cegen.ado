*! version 1.4.0 20261127 github.com/mdroste/stata-ctools
*! original: Stata's egen (egen.ado 3.4.2, April 13, 2020)
*! original: see help cegen for details

* Value label of group(varlist), label: the text of each group, from one
* observation of each group (reps, in group order). A string variable gives
* its value; a numeric one the text of its value label where that defines
* one, else the number ("%18.0g" for integers, "%9.0g" or "%9,0g"
* otherwise). The parts are joined by blanks. truncate (missing: none) cuts
* each part, by display columns when unicode is set and by characters
* otherwise; subnull writes a NUL byte as \0 instead of ending the text there.
capture mata: mata drop _cegen_group_label()
capture mata: mata drop _cegen_label_cut()
mata:
string colvector _cegen_label_cut(string colvector x, real scalar truncate, real scalar unicode)
{
    if (truncate >= .) return(x)
    if (unicode) return(ustrtrim(udsubstr(x, 1, truncate)))
    return(ustrtrim(usubstr(x, 1, truncate)))
}

void _cegen_group_label(string scalar varlist, string scalar reps, string scalar lname,
                        real scalar truncate, real scalar unicode, real scalar subnull)
{
    real rowvector v
    real colvector obs, s, nul
    real scalar i, j
    string scalar fmt, vl
    string colvector x, y

    v = st_varindex(tokens(varlist))
    obs = strtoreal(tokens(reps))'
    if (rows(obs) == 0) return
    fmt = (st_global("c(dp)") == "comma" ? "%9,0g" : "%9.0g")
    for (j = 1; j <= cols(v); j++) {
        if (st_isnumvar(v[j])) {
            s = st_data(obs, v[j])
            x = J(rows(s), 1, "")
            vl = st_varvaluelabel(v[j])
            if (vl != "") {
                if (st_vlexists(vl)) x = _cegen_label_cut(st_vlmap(vl, s), truncate, unicode)
            }
            for (i = 1; i <= rows(s); i++) {
                if (x[i] == "") x[i] = strofreal(s[i], (s[i] != trunc(s[i]) ? fmt : "%18.0g"))
            }
        }
        else {
            x = st_sdata(obs, v[j])
            if (subnull) x = usubinstr(x, char(0), "\0", .)
            else {
                nul = ustrpos(x, char(0))
                for (i = 1; i <= rows(x); i++) {
                    if (nul[i] > 0) x[i] = ustrtrim(usubstr(x[i], 1, nul[i] - 1))
                }
            }
            x = _cegen_label_cut(x, truncate, unicode)
        }
        y = (j == 1 ? x : y :+ " " :+ x)
    }
    st_vlmodify(lname, 1::rows(y), y)
}
end

program define cegen, byable(onecall)
    version 14.1
    local cvers = _caller()

    * Check observation limit (Stata plugin API limitation)
    if _N > 2147483647 {
        di as error "ctools does not support datasets exceeding 2^31 (2.147 billion) observations"
        di as error "This is a limitation of Stata's plugin API"
        exit 920
    }

    local byprefix ""
    if _by() local byprefix "`_byvars'"
    local nbyprefix : word count `byprefix'

    macro drop CTOOLS_cegen_plugin
    capture noisily _cegen_impl `cvers' `nbyprefix' `byprefix' `0'
    local rc = _rc
    if `rc' & "$CTOOLS_cegen_plugin" != "" {
        * A result the plugin still holds is not going to be stored
        capture program ctools_plugin, plugin using("$CTOOLS_cegen_plugin")
        capture plugin call ctools_plugin, "cegen discard"
    }
    macro drop CTOOLS_cegen_plugin
    exit `rc'
end

program define _cegen_impl
    version 14.1
    gettoken cvers 0 : 0
    gettoken nbyprefix 0 : 0
    local byprefix ""
    forvalues i = 1/`nbyprefix' {
        gettoken v 0 : 0
        local byprefix `byprefix' `v'
    }

    * ------------------------------------------------------------------
    * egen's grammar: [type] newvar = fcn(arguments) [if] [in] [, options]
    * ------------------------------------------------------------------
    gettoken type 0 : 0, parse(" =(")
    gettoken name 0 : 0, parse(" =(")
    local typed 1
    if `"`name'"' == "=" {
        local name `"`type'"'
        local type "`c(type)'"
        local typed 0
    }
    else {
        gettoken eqsign 0 : 0, parse(" =(")
        if `"`eqsign'"' != "=" error 198
    }
    confirm new variable `name'
    gettoken fcn 0 : 0, parse(" =(")
    gettoken args 0 : 0, parse(" ,") match(par)
    if `"`par'"' != "(" exit 198

    * Names egen kept from before Stata 9
    local fcn0 "`fcn'"
    local old   sum   any      eqany    neqany   rfirst   rlast   rmax   rmean   rmin   rmiss   robs       rsd   rsum
    local new   total anyvalue anymatch anycount rowfirst rowlast rowmax rowmean rowmin rowmiss rownonmiss rowsd rowtotal
    local pos : list posof "`fcn'" in old
    if `pos' local fcn : word `pos' of `new'

    local f_stat count mean total sd min max median pctile iqr mad mdev kurt skew std pc rank
    local f_row rowtotal rowmean rowsd rowmin rowmax rowfirst rowlast rowmiss rownonmiss ///
        rowmedian rowpctile anycount anymatch anyvalue diff
    local f_other mode seq group tag cut fill ma concat ends
    local f_all `f_stat' `f_row' `f_other'
    local byable `f_stat' mode seq

    local numtypes byte int long float double
    local isnumtype : list type in numtypes
    local isstrtype = "`type'" == "strL"
    if regexm("`type'", "^str([1-9][0-9]*)$") local isstrtype = real(regexs(1)) <= c(maxstrvarlen)
    if !`isnumtype' & !`isstrtype' {
        * Not a storage type: egen hands it to its function, whose parsing
        * decides (most stop with r(103), some ignore the type)
        _cegen_egen `cvers' "`byprefix'" `type' `name' = `fcn0'(`args') `0'
        exit
    }

    * egen parses [if] [in] [, *] before its function sees anything, and the
    * by prefix refuses an in range there (functions egen runs itself get
    * this from egen)
    if `: list fcn in f_all' {
        local rest `"`0'"'
        syntax [if] [in] [, *]
        if "`byprefix'" != "" & `"`in'"' != "" {
            di as error "'in' may not be combined with 'by'"
            exit 190
        }
        local 0 `"`rest'"'
    }

    syntax [if/] [in] [, BY(varlist) THReads(integer 0) Verbose *]
    * Later syntax statements reset these macros: keep the command's own
    local ifexp `"`if'"'
    local inrange `"`in'"'
    local byvars `by'
    local fopts `"`options'"'
    local nthreads `threads'
    local verb `verbose'
    if `nthreads' < 0 {
        di as error "cegen: threads() must be a positive integer"
        exit 198
    }
    if "`byprefix'" != "" {
        if "`byvars'" != "" {
            di as error "cegen: specify by() or the by prefix, not both"
            exit 198
        }
        local byvars `byprefix'
    }
    local typeword = cond(`typed', "`type'", "")
    local ifword ""
    if `"`ifexp'"' != "" local ifword `"if `ifexp'"'

    * A function cegen does not compute goes to egen itself: mtr(), whose tax
    * tables ctools does not copy, and community-contributed functions
    if !`: list fcn in f_all' {
        capture findfile _g`fcn'.ado
        if _rc {
            di as error "unknown {bf:egen} function {bf:`fcn'()}"
            exit 133
        }
        _cegen_egen `cvers' "`byprefix'" `typeword' `name' = `fcn0'(`args') `ifword' `inrange', ///
            `=cond("`byprefix'" == "" & "`byvars'" != "", "by(`byvars')", "")' `fopts'
        exit
    }
    * With no arguments egen's functions stop with different codes, and a
    * few create a variable: egen decides (seq() and group() take none)
    if trim(`"`args'"') == "" & !inlist("`fcn'", "seq", "group") {
        _cegen_egen `cvers' "`byprefix'" `typeword' `name' = `fcn0'(`args') `ifword' `inrange', ///
            `=cond("`byprefix'" == "" & "`byvars'" != "", "by(`byvars')", "")' `fopts'
        exit
    }
    if "`byvars'" != "" & !`: list fcn in byable' {
        di as error "{bf:cegen} ... {bf:`fcn0'()} may not be combined with {bf:by}"
        exit 190
    }

    if "`verb'" != "" quietly _ctools_clock __t_start

    * What each branch below decides
    local src ""          // source variables, in plugin order
    local fnopts ""       // function options for the plugin
    local type0 "`type'"  // the type egen starts the variable with
    local rule "asis"     // asis, promote (as replace does), auto (as compress does), str
    local varlabel ""
    local varnote ""
    local setlabel 0
    local cegen_list ""
    local cegen_punct ""
    local needs_egen 0    // inputs only egen can handle
    local keepin 0        // the in range alone selects the sample (ma)

    * ------------------------------------------------------------------
    * Statistics of an expression, by groups
    * ------------------------------------------------------------------
    if `: list fcn in f_stat' {
        local 0 `", `fopts'"'
        if "`fcn'" == "total" {
            syntax [, Missing]
        }
        else if inlist("`fcn'", "min", "max") {
            syntax [, MISSing]
        }
        else if "`fcn'" == "pctile" {
            syntax [, P(real 50) AUTOTYPE]
        }
        else if inlist("`fcn'", "median", "iqr") {
            syntax [, AUTOTYPE]
        }
        else if "`fcn'" == "std" {
            syntax [, Mean(real 0) Sd(numlist max=1 >0) STd(numlist max=1 >0)]
        }
        else if "`fcn'" == "pc" {
            syntax [, PROP]
        }
        else if "`fcn'" == "rank" {
            syntax [, Field Track Unique]
        }
        else _cegen_nooptions `fopts'

        if "`fcn'" == "pctile" {
            if `p' <= 0 | `p' >= 100 {
                di as error "{bf:p(`p')} must be between 0 and 100"
                exit 198
            }
            local fnopts "p=`p'"
        }
        if "`fcn'" == "std" {
            if "`sd'" != "" & "`std'" != "" {
                di as error "only one of {bf:sd()} or {bf:std()} is allowed"
                exit 198
            }
            if "`sd'" == "" local sd = cond("`std'" != "", "`std'", "1")
            local fnopts "smean=`mean' ssd=`sd'"
        }
        if "`fcn'" == "rank" {
            if ("`field'" != "") + ("`track'" != "") + ("`unique'" != "") > 1 {
                di as error "options {bf:field}, {bf:track} and {bf:unique} are mutually exclusive"
                exit 198
            }
            local fnopts "`field' `track' `unique'"
        }
        if "`missing'" != "" local fnopts "`fnopts' missing"
        if "`prop'" != "" local fnopts "`fnopts' prop"
        if "`autotype'" != "" {
            local rule "auto"
            local type0 "double"
        }

        * The source: a variable as it is, an expression through generate
        local e `"`args'"'
        forvalues i = 1/10 {
            gettoken inner rest : e, match(paren)
            if "`paren'" == "(" & trim(`"`rest'"') == "" local e `"`inner'"'
            else continue, break
        }
        capture unab srcvar : `e', max(1)
        if _rc == 0 {
            * egen generates a double from the expression
            capture confirm numeric variable `srcvar'
            if _rc & "`fcn'" != "count" {
                di as error "type mismatch"
                exit 109
            }
            local src `srcvar'
        }
        else {
            tempvar x
            _cegen_bysensitive `"`e'"'
            if "`fcn'" == "count" local gen `"byte `x' = 1 if !missing(`e')"'
            else local gen `"double `x' = (`e')"'
            if "`byvars'" != "" & `bysensitive' _cegen_byeval `byvars' : `gen'
            else quietly generate `gen'
            local src `x'
        }

        if "`fcn'" == "std" {
            local varlabel = `"Standardized values of "' + ustrtrim(`"`args'"')
            local setlabel 1
        }
        if "`fcn'" == "rank" {
            * egen's label holds the expression as its function received it,
            * blanks included: one before an if, two before an in, three
            * before the options, five when nothing follows
            local kind "`field'`track'`unique'"
            if "`kind'" != "" local kind "`kind' "
            local pad 5
            if `"`ifexp'"' != "" local pad 1
            else if `"`inrange'"' != "" local pad 2
            else if `"`byvars'`fopts'"' != "" local pad 3
            local varlabel = `"`kind'rank of (`args')"' + (`pad' + 1) * " "
            if "`byvars'" != "" local varlabel = `"`varlabel'"' + "by `byvars'"
            local setlabel 1
        }
    }

    * ------------------------------------------------------------------
    * mode(varname) and seq()
    * ------------------------------------------------------------------
    else if "`fcn'" == "mode" {
        local 0 `", `fopts'"'
        syntax [, MISSing MINmode MAXmode Nummode(numlist max=1 int >=1)]
        if ("`minmode'" != "") + ("`maxmode'" != "") + ("`nummode'" != "") > 1 {
            di as error "options {bf:minmode}, {bf:maxmode} and {bf:nummode()} are mutually exclusive"
            exit 198
        }
        unab src : `args', min(1) max(1)
        local type0 : type `src'
        local rule "asis"
        capture confirm string variable `src'
        if _rc == 0 {
            local rule "str"
            local type0 "str"
        }
        local fnopts "`missing' `minmode' `maxmode'"
        if "`nummode'" != "" local fnopts "`fnopts' nummode=`nummode'"
    }
    else if "`fcn'" == "seq" {
        if trim(`"`args'"') != "" error 198
        local 0 `", `fopts'"'
        syntax [, Block(integer 1) From(integer 1) To(string)]
        if `block' < 1 {
            di as error "{bf:block(`block')} must be at least 1"
            exit 498
        }
        if "`to'" == "" local to = _N
        else confirm integer number `to'
        local fnopts "from=`from' to=`to' block=`block'"
        * egen starts from a byte and lets replace widen it
        local type0 "byte"
        local rule "promote"
    }

    * ------------------------------------------------------------------
    * Functions across variables
    * ------------------------------------------------------------------
    else if `: list fcn in f_row' {
        if inlist(trim(`"`args'"'), "_all", "*") local args "_all"
        unab src : `args', min(1)
        local nsrcvars : word count `src'
        local 0 `", `fopts'"'
        if "`fcn'" == "rowtotal" {
            syntax [, Missing]
        }
        else if "`fcn'" == "rownonmiss" {
            syntax [, Strok]
        }
        else if "`fcn'" == "rowpctile" {
            syntax [, p(real 50)]
        }
        else if inlist("`fcn'", "anycount", "anymatch") {
            syntax , Values(numlist integer)
        }
        else if "`fcn'" == "anyvalue" {
            syntax , Values(string)
        }
        else _cegen_nooptions `fopts'

        local nstring 0
        foreach v of local src {
            capture confirm string variable `v'
            if _rc == 0 local ++nstring
        }
        if "`fcn'" == "rowsd" & `nsrcvars' < 2 {
            di as error "too few variables specified"
            exit 102
        }
        if "`fcn'" == "diff" & `nsrcvars' < 2 {
            di as error "too few variables specified"
            exit 102
        }
        if "`fcn'" == "anyvalue" & `nsrcvars' > 1 {
            di as error "too many variables specified"
            exit 103
        }
        if `nstring' {
            if inlist("`fcn'", "rowfirst", "rowlast") {
                * a string result, which egen builds in the type asked for
                local needs_egen 1
            }
            else if inlist("`fcn'", "rowmedian", "rowpctile") {
                foreach v of local src {
                    capture confirm numeric variable `v'
                    if _rc {
                        di as error "string variable `v' not allowed"
                        exit 109
                    }
                }
            }
            else if "`fcn'" == "rowmiss" | ("`fcn'" == "rownonmiss" & "`strok'" != "") {
            }
            else if "`fcn'" == "diff" & `nstring' == `nsrcvars' {
            }
            else {
                di as error "type mismatch"
                exit 109
            }
        }

        if "`fcn'" == "rowpctile" {
            if `p' <= 0 | `p' >= 100 {
                di as error "p(`p') must be between 0 and 100"
                exit 198
            }
            local fnopts "p=`p'"
        }
        if "`fcn'" == "rowmedian" local fnopts "p=50"
        if "`missing'" != "" local fnopts "`fnopts' missing"
        if "`strok'" != "" local fnopts "`fnopts' strok"

        if inlist("`fcn'", "anycount", "anymatch", "anyvalue") {
            local shown `"`values'"'
            if "`fcn'" == "anyvalue" {
                numlist `"`values'"', integer
                local cegen_list "`r(numlist)'"
                local label_text `"`src' if `shown'"'
                if length(`"`label_text'"') > 80 {
                    local varnote `"`label_text'"'
                    local varlabel "`src': see notes"
                }
                else local varlabel `"`label_text'"'
            }
            else {
                local cegen_list "`values'"
                local label_text "`src' == `values'"
                if length("`src'") >= 69 {
                    local varnote `"`label_text'"'
                    local varlabel "see notes"
                }
                else if length(`"`label_text'"') > 80 {
                    local varnote `"`label_text'"'
                    local varlabel "`src': see notes"
                }
                else local varlabel `"`label_text'"'
            }
            local setlabel 1
        }
        if "`fcn'" == "diff" {
            local varlabel "diff `src'"
            local setlabel 1
        }

        * egen builds these with generate and then replace, which widens an
        * integer type that cannot hold a result; the others assign once
        if inlist("`fcn'", "rowtotal", "rowmiss", "rownonmiss", "rowfirst", "rowlast") local rule "promote"
        if inlist("`fcn'", "rowmax", "rowmin", "anyvalue") local rule "promote"
        if "`fcn'" == "anycount" {
            local type0 "byte"
            local rule "promote"
        }
        if "`fcn'" == "anymatch" local type0 "byte"
    }

    * ------------------------------------------------------------------
    * group() and tag()
    * ------------------------------------------------------------------
    else if "`fcn'" == "group" {
        local src ""
        if trim(`"`args'"') != "" unab src : `args'
        local 0 `", `fopts'"'
        capture syntax [, Missing Label(string) LName(name) Truncate(numlist max=1 integer >=1) AUTOTYPE]
        local label1 ""
        if _rc {
            syntax [, Missing Label LName(name) Truncate(numlist max=1 integer >=1) AUTOTYPE]
            local label1 "`label'"
            local label ""
        }
        local replace ""
        if "`lname'" != "" {
            if `"`label'"' != "" {
                di as error "option {bf:lname()} may not be combined with {bf:label()}"
                exit 198
            }
            local label1 "label"
            local replace "replace"
        }
        if "`truncate'" != "" & "`label1'" == "" & `"`label'"' == "" {
            di as error "option {bf:truncate()} requires the {bf:label} option"
            exit 198
        }
        local wantlabel = "`label1'" != "" | `"`label'"' != ""
        if `"`label'"' != "" {
            local keep_truncate "`truncate'"
            local 0 `"`label'"'
            syntax [name(name=lname)] [, Truncate(numlist max=1 integer >=1) REPLACE]
            if "`truncate'" == "" local truncate "`keep_truncate'"
        }
        if `wantlabel' {
            if "`lname'" == "" local lname "`name'"
            if "`replace'" == "" {
                capture label list `lname'
                if _rc == 0 {
                    di as error "label {bf:`lname'} already defined"
                    di as error "    Specify {bf:label(`lname', replace)} to replace it, or name another label."
                    exit 110
                }
            }
        }
        local subnull 0
        local unicode 1
        if `cvers' <= 16 {
            local subnull 1
            local unicode 0
            if "`truncate'" == "" local truncate 32
        }
        if "`truncate'" == "" local truncate .
        local fnopts "`missing'"
        if `wantlabel' & "`src'" != "" local fnopts "`fnopts' label"
        if "`autotype'" != "" {
            local type0 "byte"
            local rule "promote"
        }
        else if !inlist("`type'", "float", "double") & `cvers' <= 16.1 local rule "promote"
        if "`src'" != "" {
            if udstrlen("group(`src')") > 80 {
                local varnote "group(`src')"
                local varlabel "see notes"
            }
            else local varlabel "group(`src')"
            local setlabel 1
        }
    }
    else if "`fcn'" == "tag" {
        unab src : `args', min(1)
        local 0 `", `fopts'"'
        syntax [, Missing]
        local fnopts "`missing'"
        local type0 "byte"
        if length("tag(`src')") > 80 local varnote "tag(`src')"
        local varlabel "tag(`src')"
        local setlabel 1
    }

    * ------------------------------------------------------------------
    * cut(), fill() and ma()
    * ------------------------------------------------------------------
    else if "`fcn'" == "cut" {
        capture unab src : `args', min(1) max(1)
        if _rc {
            di as error `"cut(`args') invalid, variable name expected"'
            exit 198
        }
        confirm numeric variable `src'
        local 0 `", `fopts'"'
        syntax [, AT(numlist ascending) Group(real 0) ICodes Label]
        if "`at'" == "" & `group' == 0 {
            di as error "one of at() or group() must be specified"
            exit 198
        }
        if "`at'" != "" & `group' != 0 {
            di as error "cannot specify both at() and group() option"
            exit 198
        }
        if `group' != 0 {
            if `group' != int(`group') | `group' < 2 {
                di as error "group() must be an integer greater than or equal to 2"
                exit 198
            }
            * pctile's limit: _N + 1, whatever the size of the sample
            if `group' > _N + 1 {
                di as error "group() must be less than or equal to the number of observations plus one"
                exit 198
            }
            local icodes "icodes"
            local fnopts "nq=`group'"
        }
        else local cegen_list "`at'"
        if "`label'" != "" local icodes "icodes"
        local fnopts "`fnopts' `icodes'"
        * egen compares a working copy of the variable in the default type
        if "`c(type)'" == "float" local fnopts "`fnopts' xfloat"
        local rule "promote"
    }
    else if "`fcn'" == "fill" {
        if `"`ifexp'`inrange'"' != "" | `"`fopts'"' != "" error 198
        numlist `"`args'"', min(2)
        local cegen_list "`r(numlist)'"
        if `: word count `cegen_list'' > _N {
            di as error "fill() lists more values than there are observations"
            exit 198
        }
        if "`c(type)'" == "float" local fnopts "xfloat"
        local rule "promote"
    }
    else if "`fcn'" == "ma" {
        if `"`ifexp'"' != "" {
            di as error "if not allowed"
            exit 101
        }
        local 0 `", `fopts'"'
        syntax [, T(integer 3) noMiss]
        if `t' < 3 | mod(`t', 2) == 0 {
            di as error "t() must be odd and >= 3"
            exit 198
        }
        local fnopts "window=`t'"
        if "`miss'" != "" local fnopts "`fnopts' nomiss"
        capture unab src : `args', max(1)
        if _rc == 0 confirm numeric variable `src'
        else {
            tempvar x
            quietly generate double `x' = (`args') `inrange'
            local src `x'
        }
        local rule "promote"
        local keepin 1
    }

    * ------------------------------------------------------------------
    * concat() and ends()
    * ------------------------------------------------------------------
    else if "`fcn'" == "concat" {
        if inlist(trim(`"`args'"'), "_all", "*") local args "_all"
        unab src : `args', min(1)
        local 0 `", `fopts'"'
        syntax [, Punct(string) Decode MAXLength(string) Format(string)]
        if "`maxlength'" != "" {
            capture confirm integer number `maxlength'
            if _rc local maxlength 0
            if `maxlength' < 1 | `maxlength' > c(maxvlabellen) {
                di as error "invalid maxlength( )"
                exit 198
            }
            local maxlength "maxlength(`maxlength')"
        }
        if "`format'" != "" {
            tempvar probe
            quietly generate byte `probe' = 1
            capture format `probe' `format'
            if _rc == 120 {
                di as error "Invalid format"
                exit 120
            }
            drop `probe'
        }
        local cegen_punct `"`punct'"'
        local rule "str"
        local type0 "str"
        * Only Stata formats a number in a display format or decodes it
        if "`format'" != "" | "`decode'" != "" {
            local converted ""
            foreach v of local src {
                capture confirm string variable `v'
                if _rc == 0 {
                    local converted `converted' `v'
                    continue
                }
                tempvar text
                _cegen_tostring `v' `text' `"`format'"' "`decode'" "`maxlength'"
                local converted `converted' `text'
            }
            local src `converted'
        }
    }
    else if "`fcn'" == "ends" {
        unab src : `args', min(1) max(1)
        capture confirm string variable `src'
        if _rc {
            di as error "type mismatch"
            exit 109
        }
        local 0 `", `fopts'"'
        syntax [, Punct(string) Head Last Tail TRim]
        if ("`head'" != "") + ("`tail'" != "") + ("`last'" != "") > 1 {
            di as error "you can only specify one option, head, tail or last"
            exit 198
        }
        if `"`punct'"' == "" local punct " "
        local cegen_punct `"`punct'"'
        local fnopts "`head' `last' `tail' `trim'"
        local rule "str"
        local type0 "str"
    }

    * mode(), seq(), tag(), anycount() and anymatch() choose their own type;
    * what egen runs itself gets its own check there
    if "`rule'" != "str" & !`isnumtype' & !`needs_egen' & ///
        !inlist("`fcn'", "mode", "seq", "tag", "anycount", "anymatch") {
        di as error "type mismatch"
        exit 109
    }

    * ------------------------------------------------------------------
    * Empty dataset: nothing to compute
    * ------------------------------------------------------------------
    if _N == 0 & !`needs_egen' {
        if "`rule'" == "str" quietly generate str1 `name' = ""
        else quietly generate `type0' `name' = .
        if `setlabel' label variable `name' `"`varlabel'"'
        exit
    }

    * ------------------------------------------------------------------
    * The sample
    * ------------------------------------------------------------------
    local marker ""
    local select "skipif"
    if `"`ifexp'"' != "" {
        local select "touse"
        * an if that names a numeric variable is itself the marker
        local bare = trim(`"`ifexp'"')
        capture confirm numeric variable `bare', exact
        if _rc == 0 & `: word count `bare'' == 1 local marker `bare'
        else {
            tempvar touse
            _cegen_bysensitive `"`ifexp'"'
            if "`byvars'" != "" & `"`inrange'"' == "" & `bysensitive' {
                _cegen_byeval `byvars' : byte `touse' = ((`ifexp') != 0)
            }
            else mark `touse' if `ifexp' `inrange'
            local marker `touse'
        }
    }

    * ------------------------------------------------------------------
    * Compute
    * ------------------------------------------------------------------
    if `needs_egen' {
        _cegen_egen `cvers' "`byprefix'" `typeword' `name' = `fcn0'(`args') `ifword' `inrange', ///
            `=cond("`byprefix'" == "" & "`byvars'" != "", "by(`byvars')", "")' `fopts'
        exit
    }

    _ctools_load
    * Stata scopes plugin registrations to the calling ado program.
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    capture confirm number 0
    global CTOOLS_cegen_plugin "`__ctools_plugin'"

    local threads_opt ""
    if `nthreads' > 0 local threads_opt "threads(`nthreads')"
    local stats = cond(inlist("`rule'", "promote", "auto"), "stats", "")
    local listlen : length local cegen_list
    mata: st_local("punctlen", strofreal(strlen(st_local("cegen_punct"))))
    if "`verb'" != "" quietly _ctools_clock __t_plugin

    * The loop is left through these flags: exit does not leave a program
    * from inside a loop
    local plugin_rc 0
    local use_egen 0
    forvalues attempt = 1/2 {
        local nsrc : word count `src'
        local nby : word count `byvars'
        local plugvars `src' `byvars' `marker'
        local __ctools_strw ""
        local anystrl 0
        foreach v of local plugvars {
            if "`: type `v''" == "strL" local anystrl 1
        }
        if "`plugvars'" != "" _ctools_strw `plugvars'
        local cegen_needfmt ""
        * A strL value the plugin cannot read (binary, or longer than 2,045
        * bytes) comes back as r(109), without a message, and is egen's
        local show = cond(`anystrl', "", "noisily")
        capture `show' plugin call ctools_plugin `plugvars' `inrange', ///
            "cegen `threads_opt' compute fn=`fcn' nsrc=`nsrc' nby=`nby' type=`type0' `select' `stats' listlen=`listlen' punctlen=`punctlen' `fnopts'"
        local plugin_rc = _rc
        if `plugin_rc' == 109 & `anystrl' {
            local use_egen 1
            local plugin_rc 0
        }
        if `plugin_rc' | `use_egen' continue, break
        if "`fcn'" != "concat" | !strpos("`cegen_needfmt'", "1") continue, break
        * concat(): numbers that are not plain integers get Stata's string()
        local converted ""
        local i 0
        foreach v of local src {
            local ++i
            if substr("`cegen_needfmt'", `i', 1) == "1" {
                tempvar text
                _cegen_tostring `v' `text' "" "" ""
                local converted `converted' `text'
            }
            else local converted `converted' `v'
        }
        local src `converted'
    }
    if `use_egen' {
        _cegen_egen `cvers' "`byprefix'" `typeword' `name' = `fcn0'(`args') `ifword' `inrange', ///
            `=cond("`byprefix'" == "" & "`byvars'" != "", "by(`byvars')", "")' `fopts'
        exit
    }
    if `plugin_rc' {
        if `plugin_rc' == 109 di as error "type mismatch"
        exit `plugin_rc'
    }

    * ------------------------------------------------------------------
    * What the function reported
    * ------------------------------------------------------------------
    if "`fcn'" == "ma" & `cegen_allmiss' {
        di as error "missing values encountered"
        exit 412
    }
    if "`fcn'" == "cut" & `cegen_badgroups' {
        di as error "group() must be less than or equal to the number of observations plus one"
        exit 198
    }
    if "`fcn'" == "cut" & `cegen_allmiss' {
        * group() over a sample with no nonmissing values: pctile's error
        di as error "no observations"
        exit 2000
    }
    if "`fcn'" == "group" & `cegen_overflow' {
        di as error "too many values for {bf:label define}"
        exit 134
    }
    if "`fcn'" == "mode" & `cegen_n' > 0 {
        if `cegen_allmiss' {
            di as text "{p 0 9 2}warning: {bf:`src'} is missing in every observation, so the mode is"
            di as text "missing. Option {bf:missing} treats missing values as a category.{p_end}"
        }
        else {
            if `cegen_groupmiss' & "`byvars'" != "" {
                di as text "{p 0 9 2}warning: {bf:`src'} is missing in every observation of at least one"
                di as text "group, where the mode is missing. Option {bf:missing} treats missing"
                di as text "values as a category.{p_end}"
            }
            if "`nummode'" != "" {
                if `nummode' > `cegen_maxmodes' {
                    if "`byvars'" == "" {
                        di as error "{bf:nummode(`nummode')} is greater than the number of modes"
                        di as error "    There " plural(`cegen_maxmodes', "is", "are") " only `cegen_maxmodes' " plural(`cegen_maxmodes', "mode") "."
                    }
                    else {
                        di as error "{bf:nummode(`nummode')} is greater than the number of modes in every group"
                        di as error "    The largest number of modes in a group is `cegen_maxmodes'."
                    }
                    exit 198
                }
                if `cegen_badgroups' {
                    di as text "{p 0 9 2}warning: {bf:nummode(`nummode')} is greater than the number of modes in"
                    di as text "`cegen_badgroups' " plural(`cegen_badgroups', "group") ", where the mode is missing.{p_end}"
                }
            }
            else if "`minmode'`maxmode'" == "" & `cegen_multiple' {
                if "`byvars'" == "" {
                    di as text "{p 0 9 2}warning: there is more than one mode, so the mode is missing."
                }
                else {
                    di as text "{p 0 9 2}warning: one or more groups have more than one mode, and the mode"
                    di as text "is missing there."
                }
                di as text "Options {bf:minmode}, {bf:maxmode} and {bf:nummode()} select one of the modes.{p_end}"
            }
        }
    }

    * ------------------------------------------------------------------
    * The new variable's storage type
    * ------------------------------------------------------------------
    local setformat ""
    if "`cegen_kind'" == "str" {
        if `cegen_strlen' > 2045 {
            * a strL result: egen's
            capture plugin call ctools_plugin, "cegen discard"
            _cegen_egen `cvers' "`byprefix'" `typeword' `name' = `fcn0'(`args') `ifword' `inrange', ///
                `=cond("`byprefix'" == "" & "`byvars'" != "", "by(`byvars')", "")' `fopts'
            exit
        }
        local ftype "str`=max(1, `cegen_strlen')'"
        local stype "str"
    }
    else {
        local ftype "`type0'"
        if "`cegen_ftype'" != "" {
            * the plugin followed egen's generate and replace steps
            local ftype "`cegen_ftype'"
            if "`ftype'" != "`type0'" local setformat = cond("`type0'" == "long", "%12.0g", "%8.0g")
        }
        else if "`rule'" == "promote" & `cegen_nonmissing' > 0 & !inlist("`type0'", "float", "double") {
            if `cegen_isint' {
                local ftype "double"
                local found 0
                foreach t in byte int long {
                    if "`t'" == "`type0'" local found 1
                    if !`found' continue
                    local lo = cond("`t'" == "byte", -127, cond("`t'" == "int", -32767, -2147483647))
                    local hi = cond("`t'" == "byte", 100, cond("`t'" == "int", 32740, 2147483620))
                    if `cegen_min' >= `lo' & `cegen_max' <= `hi' {
                        local ftype "`t'"
                        continue, break
                    }
                }
            }
            else local ftype = cond("`type0'" == "long", "double", "float")
            * replace keeps the display format the variable started with
            if "`ftype'" != "`type0'" local setformat = cond("`type0'" == "long", "%12.0g", "%8.0g")
        }
        if "`rule'" == "auto" {
            * compress: the narrowest type that holds every value exactly
            local ftype "byte"
            if `cegen_nonmissing' > 0 {
                local ftype = cond(`cegen_isfloat', "float", "double")
                if `cegen_isint' {
                    if `cegen_min' >= -127 & `cegen_max' <= 100 local ftype "byte"
                    else if `cegen_min' >= -32767 & `cegen_max' <= 32740 local ftype "int"
                    else if `cegen_min' >= -2147483647 & `cegen_max' <= 2147483620 local ftype "long"
                }
            }
            local setformat "%10.0g"
        }
        local stype "`ftype'"
    }

    * ------------------------------------------------------------------
    * Create the variable without initial values and store every observation
    * ------------------------------------------------------------------
    tempvar dummy
    mata: (void) st_addvar(st_local("ftype"), st_local("dummy"), 1)
    plugin call ctools_plugin `dummy', "cegen `threads_opt' store type=`stype'"
    if "`verb'" != "" quietly _ctools_clock __t_stored

    if "`setformat'" != "" format `dummy' `setformat'
    if `setlabel' label variable `dummy' `"`varlabel'"'
    if `"`varnote'"' != "" note `dummy' : `varnote'

    if "`fcn'" == "group" & "`lname'" != "" {
        if "`cegen_reps'" != "" {
            mata: _cegen_group_label(st_local("src"), st_local("cegen_reps"), st_local("lname"), ///
                `truncate', `unicode', `subnull')
        }
        capture label values `dummy' `lname'
    }
    if "`fcn'" == "cut" & "`label'" != "" {
        * A value label named after the variable, or the next free name
        local vlbl "`name'"
        capture label list `vlbl'
        if _rc == 0 {
            local stem = usubstr("`vlbl'", 1, ustrlen("`vlbl'") - 1)
            local tails "abcdefghijklmnopqrstuvwxyz0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
            forvalues j = 1/62 {
                local vlbl = "`stem'" + substr("`tails'", `j', 1)
                capture label list `vlbl'
                if _rc continue, break
            }
        }
        local cuts "`at'"
        if "`cuts'" == "" {
            local cuts "`cegen_cuts'"
            if `"`ifexp'`inrange'"' != "" & "`cuts'" != "" {
                * egen's first cut point is the minimum over all observations
                quietly summarize `src', meanonly
                local lowest = r(min)
                if "`c(type)'" == "float" local lowest = float(`lowest')
                gettoken first cuts : cuts
                local cuts "`lowest' `cuts'"
            }
        }
        local code 0
        foreach cutp of local cuts {
            local text = string(`cutp')
            label define `vlbl' `code' "`text'-", add
            local ++code
        }
        label values `dummy' `vlbl'
    }

    rename `dummy' `name'
    * Plugin stores do not mark the data changed; generate does
    mata: st_updata(1)

    if `cegen_nmiss' > 0 {
        local shown : display %11.0fc `cegen_nmiss'
        local shown = trim("`shown'")
        di as text "(`shown' missing value" cond(`cegen_nmiss' > 1, "s", "") " generated)"
    }

    if "`verb'" != "" {
        quietly _ctools_clock __t_end
        di as text ""
        di as text "{hline 60}"
        di as text "cegen timing breakdown:"
        di as text "{hline 60}"
        di as text "  Stata pre-plugin:         " as result %8.4f `__t_plugin' - `__t_start' " sec"
        di as text "  C plugin phases:"
        di as text "    Data load:              " as result %8.4f `cegen_t_load' " sec"
        di as text "    Compute `fcn'():" _col(29) as result %8.4f `cegen_t_compute' " sec"
        di as text "    Data store:             " as result %8.4f `cegen_t_store' " sec"
        di as text "  Stata post-plugin:        " as result %8.4f `__t_end' - `__t_stored' " sec"
        di as text "{hline 60}"
        di as text "    Wall clock total:       " as result %8.4f `__t_end' - `__t_start' " sec"
        di as text "{hline 60}"
        di as text "  Sample:  " as result `cegen_n' as text " of " as result _N as text " observations"
        di as text "  Groups:  " as result `cegen_ngroups'
        di as text "  Type:    " as result "`ftype'"
        di as text "  Threads: " as result `cegen_threads'
        di as text "{hline 60}"
    }
end

* 1 in the caller's local bysensitive when an expression depends on the
* by-group it is evaluated in: _n, _N, an explicit subscript, or sum()
program define _cegen_bysensitive
    version 14.1
    args e
    local hit = regexm(`"`e'"', "(^|[^A-Za-z0-9_])(_n|_N)($|[^A-Za-z0-9_])") | ///
        strpos(`"`e'"', "[") > 0 | regexm(`"`e'"', "(^|[^A-Za-z0-9_])sum[ ]*\(")
    c_local bysensitive `hit'
end

* `by varlist: generate ...' for a statement that depends on the by-group.
* Usage: _cegen_byeval varlist : type newvar = exp
program define _cegen_byeval, sortpreserve
    version 14.1
    gettoken by 0 : 0, parse(":")
    gettoken colon 0 : 0, parse(":")
    sort `by'
    quietly by `by': generate `0'
end

* The text of a numeric variable as concat() takes it, in new variable text:
* its value label (decode) where one applies, else string() in the display
* format given
program define _cegen_tostring
    version 14.1
    args v text format decode maxlength
    local fmt ""
    if `"`format'"' != "" local fmt `", "`format'""'
    local done 0
    if "`decode'" != "" {
        capture decode `v', generate(`text') `maxlength'
        if _rc == 0 {
            quietly replace `text' = string(`v'`fmt') if `text' == ""
            local done 1
        }
        else capture drop `text'
    }
    if !`done' quietly generate `text' = string(`v'`fmt')
end

* Stop when a function that takes no options was given some
program define _cegen_nooptions
    version 14.1
    if `"`0'"' == "" exit
    gettoken bad : 0, parse(" (")
    di as error "option {bf:`bad'} not allowed"
    exit 198
end

* Run egen itself, under the caller's version and by prefix
program define _cegen_egen
    gettoken cvers 0 : 0
    gettoken byprefix 0 : 0
    if "`byprefix'" != "" {
        version `cvers', missing: by `byprefix': egen `0'
    }
    else {
        version `cvers', missing: egen `0'
    }
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
