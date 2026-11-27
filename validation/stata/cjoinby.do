* cjoinby: named cases share this file, each do-call has its own local scope.
gettoken __case 0 : 0
macro shift
if "`__case'" == "" local __case "cjoinby"

if "`__case'" == "cjoinby" {
* Reference tests compare with Stata's joinby (official).
do "validation/support/validate_setup.do"
capture confirm file "validation/support/rollback_helpers.do"
if !_rc do "validation/support/rollback_helpers.do"
else do "rollback_helpers.do"
which joinby

* _cjoinby_parity: run joinby and cjoinby on the same master and using data
* and compare return codes and, when both succeed, everything joinby leaves:
* the observations (joinby's order within ties depends on the sort seed, so
* both results are compared in the order of all variables), variable order,
* storage types, formats, variable and value labels, value-label definitions,
* notes and other characteristics, the data label, the sort marker and
* c(changed) (_ctools_rb_state). A sort marker must also hold in cjoinby's
* own order.
capture program drop _cjoinby_parity
program define _cjoinby_parity
    syntax, master(string) using(string) testname(string) [keys(string) ///
        opts(string asis) copts(string asis)]
    tempfile ref
    quietly use `"`master'"', clear
    capture joinby `keys' using `"`using'"', `opts'
    local rc_native = _rc
    if !`rc_native' {
        _ctools_rb_state
        local st_native `"`r(state)'"'
        local all
        foreach v of varlist _all {
            local all `all' `v'
        }
        sort `all'
        quietly save `"`ref'"'
    }
    quietly use `"`master'"', clear
    capture cjoinby `keys' using `"`using'"', `opts' `copts'
    local rc_ctools = _rc
    if `rc_native' != `rc_ctools' {
        test_fail "`testname'" "joinby r(`rc_native'), cjoinby r(`rc_ctools')"
        exit
    }
    if `rc_native' {
        test_pass "`testname'"
        exit
    }
    _ctools_rb_state
    local st_ctools `"`r(state)'"'
    mata: st_local("same", strofreal(st_local("st_native") == st_local("st_ctools")))
    if !`same' {
        mata: _cjoinby_first_difference(st_local("st_native"), st_local("st_ctools"))
        test_fail "`testname'" `"state differs: joinby [`macval(d_native)'] cjoinby [`macval(d_ctools)']"'
        exit
    }
    local sortedby : sortedby
    if "`sortedby'" != "" & _N > 0 {
        tempvar pos
        quietly generate long `pos' = _n
        sort `sortedby', stable
        capture assert `pos' == _n
        if _rc {
            test_fail "`testname'" "the data are not in the order of the sort marker (`sortedby')"
            exit
        }
        drop `pos'
    }
    sort `all'
    capture cf _all using `"`ref'"', all
    if _rc test_fail "`testname'" "observations differ from joinby (cf r(`=_rc'))"
    else test_pass "`testname'"
end

capture mata: mata drop _cjoinby_first_difference()
mata:
/* Where two state strings first differ, with some context, in locals
   d_native and d_ctools */
void _cjoinby_first_difference(string scalar a, string scalar b)
{
    real scalar i, n

    n = min((strlen(a), strlen(b)))
    for (i = 1; i <= n; i++) {
        if (substr(a, i, 1) != substr(b, i, 1)) break
    }
    st_local("d_native", substr(a, max((1, i - 60)), 140))
    st_local("d_ctools", substr(b, max((1, i - 60)), 140))
}
end

* _cjoinby_data <file> <n> <seed> <m|u> <keytype>: a master (m) or using (u)
* dataset with duplicated and missing keys, two variables both datasets have
* (sh, shs), variables of its own, labels, formats, notes and characteristics.
capture program drop _cjoinby_data
program define _cjoinby_data
    args file n seed role keytype
    clear
    set seed `seed'
    quietly set obs `n'
    if "`role'" == "m" gen long mid = _n
    else gen long uid = _n
    if "`keytype'" == "int" {
        gen int id = floor(runiform() * 12)
        replace id = . if runiform() < .08
        replace id = .a if runiform() < .05
    }
    else if "`keytype'" == "double" {
        gen double id = floor(runiform() * 12) / 4 + .1
        replace id = .b if runiform() < .08
        replace id = -0 if runiform() < .05
    }
    else if "`keytype'" == "sparse" {
        gen long id = 100000000 + floor(runiform() * 12) * 7919
        replace id = -5 if runiform() < .05
    }
    else if "`keytype'" == "string" {
        gen str8 id = "k" + string(floor(runiform() * 12))
        replace id = "" if runiform() < .08
        replace id = "k1 " if runiform() < .05
    }
    else if "`keytype'" == "two" {
        gen byte g = floor(runiform() * 4)
        gen int id = floor(runiform() * 5) - 2
        replace g = . if runiform() < .06
        replace id = .z if runiform() < .06
    }
    else if "`keytype'" == "mixed" {
        gen str3 g = "g" + string(floor(runiform() * 3))
        gen double id = floor(runiform() * 5) + .5
        replace g = "" if runiform() < .06
    }
    else error 198
    gen double sh = floor(runiform() * 3)
    replace sh = . if runiform() < .3
    replace sh = .c if runiform() < .1
    gen str5 shs = "s" + string(floor(runiform() * 3))
    replace shs = "" if runiform() < .3
    if "`role'" == "m" {
        gen double mx = runiform()
        gen str4 ms = "m" + string(mod(_n, 7))
        label variable sh "master sh"
        label define mlab 0 "m zero" 1 "m one"
        label values mid mlab
        note: master dataset note
        note sh: master note on sh
        char _dta[rb_char] master characteristic
        char sh[both] master sh characteristic
        label data "cjoinby master"
    }
    else {
        gen double ux = runiform()
        gen str6 us = "u" + string(mod(_n, 5))
        gen byte uk = mod(_n, 3)
        label variable sh "using sh"
        label variable ux "using ux"
        label variable id "using id"
        label define ulab 0 "u zero" 1 "u one" 2 "u two"
        label define mlab 0 "using definition of mlab"
        label values uk ulab
        label values uid mlab
        format ux %12.4f
        format sh %6.1f
        note: using dataset note
        note: master dataset note
        note sh: using note on sh
        note ux: using note on ux
        note id: using note on id
        char _dta[rb_uchar] using characteristic
        char sh[both] using sh characteristic
        char sh[usingonly] using sh characteristic 2
        char ux[own] ux characteristic
        char id[keychar] id characteristic
        label data "cjoinby using"
    }
    quietly save `"`file'"', replace
end

tempfile m u ref out

* ------------------------------------------------------------------------
print_section "Parity with joinby: key types and options"
* ------------------------------------------------------------------------
foreach keytype in int double sparse string two mixed {
    local keys id
    if inlist("`keytype'", "two", "mixed") local keys g id
    _cjoinby_data `m' 60 1101 m `keytype'
    _cjoinby_data `u' 45 2202 u `keytype'
    foreach unmatched in none both master using {
        foreach how in "" "update" "update replace" {
            _cjoinby_parity, master(`m') using(`u') keys(`keys') ///
                opts(unmatched(`unmatched') `how') ///
                testname("`keytype' keys, unmatched(`unmatched') `how'")
        }
    }
    * the default, a named _merge, nolabel and the common-variable form
    _cjoinby_parity, master(`m') using(`u') keys(`keys') testname("`keytype' keys, defaults")
    _cjoinby_parity, master(`m') using(`u') keys(`keys') opts(_merge(src)) ///
        testname("`keytype' keys, _merge(src)")
    _cjoinby_parity, master(`m') using(`u') keys(`keys') opts(nolabel unmatched(both)) ///
        testname("`keytype' keys, nolabel")
    _cjoinby_parity, master(`m') using(`u') keys(`keys') opts(nolabel update unm(b)) ///
        testname("`keytype' keys, nolabel update, abbreviated unmatched()")
}

* Common variables: joinby joins on every variable the datasets share
quietly {
    use `u', clear
    drop sh shs
    save `ref', replace
}
noisily _cjoinby_parity, master(`m') using(`ref') opts(unmatched(both)) ///
    testname("common variables as keys")
quietly {
    use `u', clear
    keep uid ux
    rename uid other
    save `ref', replace
}
_cjoinby_parity, master(`m') using(`ref') testname("no common variables (r(198))")

* ------------------------------------------------------------------------
print_section "Parity with joinby: sorted masters"
* ------------------------------------------------------------------------
* joinby is a sortpreserve program: the result of a sorted master is sorted
* the same way, with the using-only observations placed by their values of
* the sort variables.
_cjoinby_data `m' 70 3303 m two
_cjoinby_data `u' 50 4404 u two
foreach sortvars in "g id" "id g" "g" "id" "g id mid" "mx" "ms" "g ms" "ms g" "sh" "shs id" "g sh mid" "id shs" {
    quietly use `m', clear
    sort `sortvars'
    quietly save `ref', replace
    foreach unmatched in none both using master {
        foreach how in "" "update replace" {
            _cjoinby_parity, master(`ref') using(`u') keys(g id) ///
                opts(unmatched(`unmatched') `how') ///
                testname("master sorted by `sortvars', unmatched(`unmatched') `how'")
        }
    }
}
_cjoinby_data `m' 70 5505 m string
_cjoinby_data `u' 50 6606 u string
foreach sortvars in "id" "id mid" "ms id" "shs" {
    quietly use `m', clear
    sort `sortvars'
    quietly save `ref', replace
    foreach unmatched in none both {
        _cjoinby_parity, master(`ref') using(`u') keys(id) opts(unmatched(`unmatched')) ///
            testname("string key, master sorted by `sortvars', unmatched(`unmatched')")
    }
}

* ------------------------------------------------------------------------
print_section "Order of the observations"
* ------------------------------------------------------------------------
* cjoinby's own contract: each master observation is followed by its
* pairings in the using data's order; using-only observations follow, by key.
_cjoinby_data `m' 200 7707 m int
_cjoinby_data `u' 150 8808 u int
quietly {
    use `m', clear
    cjoinby id using `u', unmatched(both)
    gen long pos = _n
    * master observations in their order, using-only ones after them
    gen byte tail = _merge == 2
    sort tail mid uid, stable
    capture assert pos == _n if !tail
    local rc1 = _rc
    sort pos
    capture assert tail[_n - 1] <= tail if _n > 1
    local rc2 = _rc
    * using-only observations in key order, then in the using data's order
    capture assert (id[_n - 1] < id | (id[_n - 1] == id & uid[_n - 1] < uid)) ///
        if tail & tail[_n - 1] & _n > 1
    local rc3 = _rc
}
if !`rc1' & !`rc2' & !`rc3' test_pass "master order, pairings in using order, using-only by key"
else test_fail "master order, pairings in using order, using-only by key" "r(`rc1') r(`rc2') r(`rc3')"

* Repeated runs and thread counts give the same observations in the same order
quietly {
    use `m', clear
    cjoinby id using `u', unmatched(both) update threads(1)
    save `ref', replace
    local same 1
    foreach t in 2 4 8 {
        use `m', clear
        cjoinby id using `u', unmatched(both) update threads(`t')
        capture cf _all using `ref', all
        if _rc local same 0
    }
}
if `same' test_pass "threads(1), (2), (4) and (8) agree row for row"
else test_fail "threads(1), (2), (4) and (8) agree row for row" "results differ"

* ------------------------------------------------------------------------
print_section "Storage types"
* ------------------------------------------------------------------------
foreach mt in byte int long float double {
    foreach ut in byte int long float double {
        quietly {
            clear
            set obs 6
            gen `ut' id = mod(_n, 4)
            gen `ut' sh = _n
            gen byte un = 1
            save `u', replace
            clear
            set obs 5
            gen `mt' id = _n
            gen `mt' sh = cond(mod(_n, 2), ., _n * 10)
            format id %5.0f
            save `m', replace
        }
        foreach how in "" "update" {
            _cjoinby_parity, master(`m') using(`u') keys(id) opts(unmatched(both) `how') ///
                testname("master `mt', using `ut' `how'")
        }
    }
}
foreach pair in "str3 str9" "str9 str3" "str5 str5" "str1 str40" {
    gettoken mt ut : pair
    local ut = strtrim("`ut'")
    quietly {
        clear
        set obs 6
        gen `ut' id = substr("abcdefgh", 1, mod(_n, 3) + 1)
        gen `ut' sh = "u"
        gen byte un = 1
        save `u', replace
        clear
        set obs 5
        gen `mt' id = substr("abcdefgh", 1, mod(_n, 3) + 1)
        gen `mt' sh = cond(mod(_n, 2), "", "m")
        format id %-12s
        save `m', replace
    }
    foreach how in "" "update" "update replace" {
        _cjoinby_parity, master(`m') using(`u') keys(id) opts(unmatched(both) `how') ///
            testname("master `mt', using `ut' `how'")
    }
}
* A shared variable that is a string in one dataset and a number in the other
quietly {
    clear
    set obs 4
    gen int id = _n
    gen str3 sh = "abc"
    gen double sh2 = 7
    save `u', replace
    clear
    set obs 4
    gen int id = _n + 1
    gen double sh = cond(_n == 1, ., 5)
    gen str2 sh2 = cond(_n == 2, "", "zz")
    save `m', replace
}
foreach how in "" "update" "update replace" {
    _cjoinby_parity, master(`m') using(`u') keys(id) opts(unmatched(both) `how') ///
        testname("string/numeric shared variable `how'")
}
* float and double keys compare as doubles
quietly {
    clear
    set obs 4
    gen double id = _n + .1
    gen u = _n
    save `u', replace
    clear
    set obs 4
    gen float id = _n + .1
    replace id = 2.5 in 4
    save `m', replace
    use `u', clear
    replace id = 2.5 in 1
    save `u', replace
}
_cjoinby_parity, master(`m') using(`u') keys(id) opts(unmatched(both)) ///
    testname("float master key against double using key")

* ------------------------------------------------------------------------
print_section "Many variables"
* ------------------------------------------------------------------------
* Hundreds of shared and new variables of mixed types: the per-variable
* bookkeeping (storage types, payload positions, metadata) must scale and
* keep every variable in its place.
quietly {
    clear
    set seed 606
    set obs 90
    gen long uid = _n
    gen int id = floor(runiform() * 30)
    forvalues j = 1/150 {
        if mod(`j', 3) == 0 gen str6 sh`j' = cond(runiform() < .4, "", "u`j'")
        else if mod(`j', 3) == 1 gen double sh`j' = cond(runiform() < .4, ., `j' + .5)
        else gen long sh`j' = cond(runiform() < .4, .a, `j' * 1000)
        label variable sh`j' "using sh`j'"
    }
    forvalues j = 1/150 {
        if mod(`j', 4) == 0 gen str`=1 + mod(`j', 9)' u`j' = substr("abcdefghij", 1, mod(_n, 1 + mod(`j', 9)))
        else gen float u`j' = `j' + _n
        note u`j': using note `j'
    }
    save `u', replace
    clear
    set obs 120
    gen long mid = _n
    gen int id = floor(runiform() * 36)
    forvalues j = 1/150 {
        if mod(`j', 3) == 0 gen str3 sh`j' = cond(runiform() < .4, "", "m")
        else if mod(`j', 3) == 1 gen float sh`j' = cond(runiform() < .4, ., `j')
        else gen byte sh`j' = cond(runiform() < .4, ., mod(`j', 100))
    }
    forvalues j = 1/150 {
        gen int m`j' = `j' - _n
    }
    save `m', replace
}
foreach how in "" "update" "update replace" {
    _cjoinby_parity, master(`m') using(`u') keys(id) opts(unmatched(both) `how') ///
        testname("300 variables on each side `how'")
}

* ------------------------------------------------------------------------
print_section "strL variables"
* ------------------------------------------------------------------------
quietly {
    clear
    set seed 911
    set obs 40
    gen long uid = _n
    gen int id = floor(runiform() * 9)
    gen strL ubig = "using strL " + string(_n) + " " + 3000 * "x"
    replace ubig = "" if mod(_n, 7) == 0
    gen strL both = cond(mod(_n, 3) == 0, "", "u" + string(_n))
    gen str4 mixed = cond(mod(_n, 4) == 0, "", "u" + string(mod(_n, 9)))
    gen strL mixed2 = cond(mod(_n, 5) == 0, "", "U" + string(_n))
    label variable ubig "using strL"
    note ubig: note on a using strL
    save `u', replace
    clear
    set obs 55
    gen long mid = _n
    gen int id = floor(runiform() * 11)
    gen strL mbig = "master strL " + string(_n)
    replace mbig = "" if mod(_n, 6) == 0
    gen strL both = cond(mod(_n, 2) == 0, "", "m" + string(_n))
    gen strL mixed = cond(mod(_n, 3) == 0, "", "m" + string(_n))
    gen str5 mixed2 = cond(mod(_n, 2) == 0, "", "m" + string(mod(_n, 9)))
    label variable mbig "master strL"
    format mbig %30s
    note mbig: note on a master strL
    char mbig[own] master strL characteristic
    save `m', replace
}
foreach unmatched in none both master using {
    foreach how in "" "update" "update replace" {
        _cjoinby_parity, master(`m') using(`u') keys(id) opts(unmatched(`unmatched') `how') ///
            testname("strL variables, unmatched(`unmatched') `how'")
    }
}
quietly {
    use `m', clear
    sort id mid
    save `ref', replace
}
_cjoinby_parity, master(`ref') using(`u') keys(id) opts(unmatched(both) update) ///
    testname("strL variables, sorted master")
quietly {
    use `m', clear
    sort mbig
    save `ref', replace
}
_cjoinby_parity, master(`ref') using(`u') keys(id) opts(unmatched(both)) ///
    testname("strL variables, master sorted by a strL")
* an m:1 join moves no master observation
quietly {
    use `u', clear
    bysort id (uid): keep if _n == 1
    save `ref', replace
}
foreach how in "" "update replace" {
    _cjoinby_parity, master(`m') using(`ref') keys(id) opts(unmatched(master) `how') ///
        testname("strL variables, no master observation moves `how'")
}
* strL keys are refused, as merge refuses them
quietly {
    clear
    set obs 3
    gen strL id = "k" + string(_n)
    gen x = 1
    save `u', replace
    clear
    set obs 3
    gen strL id = "k" + string(_n)
    save `m', replace
    clear
    set obs 3
    gen str5 id = "k" + string(_n)
    save `ref', replace
}
_cjoinby_parity, master(`m') using(`u') keys(id) testname("strL key in the master (r(106))")
_cjoinby_parity, master(`ref') using(`u') keys(id) testname("strL key in the using data (r(106))")

* ------------------------------------------------------------------------
print_section "Wide strings (column transport)"
* ------------------------------------------------------------------------
quietly {
    clear
    set seed 77
    set obs 30
    gen long uid = _n
    gen int id = floor(runiform() * 6)
    gen str2045 uwide = 2040 * "u" + string(_n)
    gen str2045 shw = cond(mod(_n, 3) == 0, "", 2000 * "s" + string(_n))
    save `u', replace
    clear
    set obs 25
    gen long mid = _n
    gen int id = floor(runiform() * 8)
    gen str2045 mwide = 2040 * "m" + string(_n)
    gen str2045 shw = cond(mod(_n, 2) == 0, "", 1990 * "t" + string(_n))
    save `m', replace
}
foreach how in "" "update" "update replace" {
    _cjoinby_parity, master(`m') using(`u') keys(id) opts(unmatched(both) `how') ///
        testname("str2045 variables `how'")
}
quietly {
    clear
    set obs 20
    gen str2045 id = 2030 * "k" + string(mod(_n, 5))
    gen long uid = _n
    save `u', replace
    clear
    set obs 20
    gen str2045 id = 2030 * "k" + string(mod(_n, 7))
    gen long mid = _n
    save `m', replace
}
_cjoinby_parity, master(`m') using(`u') keys(id) opts(unmatched(both)) testname("str2045 key")

* ------------------------------------------------------------------------
print_section "Empty data and results"
* ------------------------------------------------------------------------
_cjoinby_data `m' 30 1212 m int
_cjoinby_data `u' 20 3434 u int
quietly {
    use `m', clear
    keep if 0
    save `ref', replace
}
foreach unmatched in none both master using {
    _cjoinby_parity, master(`ref') using(`u') keys(id) opts(unmatched(`unmatched')) ///
        testname("empty master, unmatched(`unmatched')")
}
* an empty master still has a sort marker, which the result must satisfy
foreach sortvars in "id" "sh id" "ms" {
    quietly {
        use `m', clear
        keep if 0
        sort `sortvars'
        save `ref', replace
    }
    foreach how in "" "update" {
        _cjoinby_parity, master(`ref') using(`u') keys(id) opts(unmatched(using) `how') ///
            testname("empty master sorted by `sortvars', unmatched(using) `how'")
    }
}
quietly {
    use `u', clear
    replace id = id + 1000
    save `ref', replace
}
foreach unmatched in none both master using {
    foreach how in "" "update" {
        _cjoinby_parity, master(`m') using(`ref') keys(id) opts(unmatched(`unmatched') `how') ///
            testname("no key in common, unmatched(`unmatched') `how'")
    }
}
* joinby stops with an error when the using data have no observations
* (its by-group counter is never created); cjoinby joins without matches
quietly {
    use `u', clear
    keep if 0
    save `ref', replace
    use `m', clear
    cjoinby id using `ref', unmatched(master)
    local ok = _N == 30 & r(N_1) == 30 & r(N_3) == 0
    capture confirm variable ux
    local ok = `ok' & !_rc
    capture assert _merge == 1 & mid == _n
    local ok = `ok' & !_rc
    use `m', clear
    cjoinby id using `ref'
    local ok = `ok' & _N == 0
}
if `ok' test_pass "using data without observations"
else test_fail "using data without observations" "unexpected result"

* ------------------------------------------------------------------------
print_section "Errors"
* ------------------------------------------------------------------------
_cjoinby_data `m' 30 1212 m int
_cjoinby_data `u' 20 3434 u int
_cjoinby_parity, master(`m') using(`u') keys(mx) testname("key not in the using data (r(111))")
_cjoinby_parity, master(`m') using(`u') keys(nosuchvar) testname("key not in the master (r(111))")
_cjoinby_parity, master(`m') using(`u') keys(id) opts(unmatched(sometimes)) ///
    testname("invalid unmatched() (r(198))")
_cjoinby_parity, master(`m') using(`u') keys(id) opts(_merge(uid)) ///
    testname("_merge() names a using variable (r(110))")
_cjoinby_parity, master(`m') using(`u') keys(id) opts(_merge(mid) unmatched(both)) ///
    testname("_merge() names a master variable (r(110))")
_cjoinby_parity, master(`m') using(`m'_absent) keys(id) testname("using file not found (r(601))")
quietly {
    use `m', clear
    gen byte _merge = 1
    save `ref', replace
}
_cjoinby_parity, master(`ref') using(`u') keys(id) opts(unmatched(both)) ///
    testname("_merge already defined (r(110))")
_cjoinby_parity, master(`ref') using(`u') keys(id) ///
    testname("a master _merge is no obstacle without unmatched()")
quietly {
    use `m', clear
    tostring id, replace force
    save `ref', replace
}
_cjoinby_parity, master(`ref') using(`u') keys(id) testname("string key against numeric key (r(106))")
* replace needs update; a failed call leaves the master as it was (joinby
* leaves a half-merged dataset here, so only the return code is compared)
quietly use `m', clear
capture cjoinby id using `u', replace
local rc = _rc
capture cf _all using `m', all
if `rc' == 198 & !_rc test_pass "replace without update (r(198)), master unchanged"
else test_fail "replace without update (r(198)), master unchanged" "r(`rc'), cf r(`=_rc')"
* A join past Stata's limit of 2,147,483,647 observations (here 50,000 x
* 50,000 pairings of one key) stops in the plan, before any change
quietly {
    clear
    set obs 50000
    gen byte id = 1
    gen long uid = _n
    save `ref', replace
    clear
    set obs 50000
    gen byte id = 1
    gen long mid = _n
    save `out', replace
}
capture cjoinby id using `ref'
local rc = _rc
capture cf _all using `out', all
local left : all globals "CTOOLS_cjoinby*"
if `rc' == 920 & !_rc & "`left'" == "" test_pass "output past the observation limit (r(920)), master unchanged"
else test_fail "output past the observation limit (r(920)), master unchanged" "r(`rc'), cf r(`=_rc'), globals `left'"
quietly use `m', clear
capture cjoinby id using `u', snapshot fast
if _rc == 198 test_pass "snapshot and fast together refused"
else test_fail "snapshot and fast together refused" "r(`=_rc')"
capture cjoinby id using `u', threads(-1)
if _rc == 198 test_pass "threads(-1) refused"
else test_fail "threads(-1) refused" "r(`=_rc')"

* ------------------------------------------------------------------------
print_section "Value label __MERGE"
* ------------------------------------------------------------------------
* joinby defines __MERGE even without a _merge variable, replaces its own
* earlier definition and leaves a user's __MERGE alone
quietly {
    use `m', clear
    label define __MERGE 1 "mine"
    save `ref', replace
}
_cjoinby_parity, master(`ref') using(`u') keys(id) opts(unmatched(both)) ///
    testname("a user's __MERGE stays")
quietly {
    use `m', clear
    joinby id using `u', update
    drop uid ux us uk
    save `ref', replace
}
_cjoinby_parity, master(`ref') using(`u') keys(id) opts(unmatched(both)) ///
    testname("joinby's own __MERGE is replaced")

* ------------------------------------------------------------------------
print_section "Stored results, options and file names"
* ------------------------------------------------------------------------
quietly {
    use `m', clear
    joinby id using `u', unmatched(both)
    forvalues c = 1/3 {
        count if _merge == `c'
        local n`c' = r(N)
    }
    local n = _N
    use `m', clear
    cjoinby id using `u', unmatched(both)
    local ok = r(N) == `n' & r(N_1) == `n1' & r(N_2) == `n2' & r(N_3) == `n3' & ///
        "`r(keyvars)'" == "id" & r(time) == .
    save `out', replace
    use `m', clear
    cjoinby id using `u', unmatched(both) verbose
    local ok = `ok' & r(time) < . & r(N) == `n'
    capture cf _all using `out', all
    local ok = `ok' & !_rc
}
if `ok' test_pass "r() results; verbose changes nothing"
else test_fail "r() results; verbose changes nothing" "unexpected r() or data"
foreach mode in snapshot fast {
    quietly use `m', clear
    quietly cjoinby id using `u', unmatched(both) `mode'
    capture cf _all using `out', all
    if !_rc test_pass "option `mode' gives the default results"
    else test_fail "option `mode' gives the default results" "results differ"
}
* A using file name with spaces, given without its extension
local spaced "`c(tmpdir)'/cjoinby using file `c(pid)'"
quietly {
    use `u', clear
    save "`spaced'.dta", replace
    use `m', clear
}
capture cjoinby id using "`spaced'", unmatched(both)
local rc = _rc
capture cf _all using `out', all
if !`rc' & !_rc test_pass "using file name with spaces, without .dta"
else test_fail "using file name with spaces, without .dta" "r(`rc'), cf r(`=_rc')"
capture erase "`spaced'.dta"

* ------------------------------------------------------------------------
print_section "Large joins"
* ------------------------------------------------------------------------
* Enough observations for every parallel kernel: direct and hashed lookups,
* numeric and string keys, expanding and filtering joins.
foreach kind in direct hashed string composite {
    quietly {
        clear
        set seed 20261126
        set obs 60000
        gen long uid = _n
        if "`kind'" == "direct" gen long id = floor(runiform() * 30000)
        else if "`kind'" == "hashed" gen double id = floor(runiform() * 30000) * 1000003 + .5
        else if "`kind'" == "string" gen str12 id = "id" + string(floor(runiform() * 30000), "%08.0f")
        else {
            gen int id = floor(runiform() * 300)
            gen byte g = floor(runiform() * 100)
        }
        gen double ux = runiform()
        gen str7 us = "u" + string(mod(_n, 997))
        gen float sh = cond(runiform() < .2, ., runiform())
        save `u', replace
        clear
        set obs 90000
        gen long mid = _n
        if "`kind'" == "direct" gen long id = floor(runiform() * 40000)
        else if "`kind'" == "hashed" gen double id = floor(runiform() * 40000) * 1000003 + .5
        else if "`kind'" == "string" gen str12 id = "id" + string(floor(runiform() * 40000), "%08.0f")
        else {
            gen int id = floor(runiform() * 400)
            gen byte g = floor(runiform() * 100)
        }
        gen double mx = runiform()
        gen str5 ms = "m" + string(mod(_n, 101))
        gen float sh = cond(runiform() < .2, ., runiform())
        save `m', replace
    }
    local keys id
    if "`kind'" == "composite" local keys g id
    foreach spec in "unmatched(none)" "unmatched(both) update" "unmatched(using) update replace" {
        _cjoinby_parity, master(`m') using(`u') keys(`keys') opts(`spec') ///
            testname("large join, `kind' keys, `spec'")
    }
    quietly {
        use `m', clear
        sort `keys' mid
        save `ref', replace
    }
    _cjoinby_parity, master(`ref') using(`u') keys(`keys') opts(unmatched(both)) ///
        testname("large join, `kind' keys, sorted master")
}

* ------------------------------------------------------------------------
print_section "Rollback of Failed Joins"
* ------------------------------------------------------------------------
* A failed or interrupted cjoinby leaves the master as it was without a
* snapshot of the dataset (Stata 16 or later); option snapshot preserves it
* instead and option fast undoes nothing (see validation/support/rollback_helpers.do).
* The hook $CTOOLS_cjoinby_fault fails the join:
* schema   with r(1) after the storage types were widened and the new
*          variables, value labels, notes, characteristics and observations
*          were added
* store    inside the plugin after half the output rows were written, with
*          the memory error it simulates, r(920)
* post     with r(1) after the plugin wrote the whole join
* strl     with r(1) after the strL values were copied
clear
set seed 4242
set obs 300
gen long uid = _n
gen long id = floor(runiform() * 40)
gen str9 skey = "s" + string(mod(_n, 40))
gen double sh = runiform()
gen str12 shs = "using" + string(_n)
gen double ux = _n / 7
gen strL ubig = "using strL " + string(_n)
label define rb_u 1 "first" 2 "second"
label values uid rb_u
label variable sh "using sh"
label data "Using data"
note: using dataset note
note sh: using note on sh
note id: using note on id
char _dta[rb_uchar] using characteristic
char sh[rb_uchar] using sh characteristic
tempfile rb_using
save `rb_using'
clear
set obs 400
gen long mid = _n
gen byte id = floor(runiform() * 50)
gen str3 skey = "s" + string(mod(_n, 50))
gen float sh = cond(mod(_n, 3) == 0, ., _n)
gen str4 shs = cond(mod(_n, 4) == 0, "", "m" + string(mod(_n, 99)))
gen double mx = runiform() * 270
gen strL mbig = "master strL " + string(_n)
label define rb_m 1 "one"
label values id rb_m
label variable mx "master mx"
label data "cjoinby rollback master"
note: master dataset note
note sh: master note on sh
char _dta[rb_char] master characteristic
char sh[rb_char] sh characteristic
sort id mid
tempfile rb_copy
save `rb_copy'

foreach fault in schema store post strl {
    local rc = cond("`fault'" == "store", 920, 1)
    _ctools_rb_case, cmd(cjoinby id using `rb_using', unmatched(both) update replace) ///
        global(CTOOLS_cjoinby_fault) fault(`fault') rc(`rc') testname("[rollback] `fault'")
}
_ctools_rb_case, cmd(cjoinby skey using `rb_using', unmatched(using)) ///
    global(CTOOLS_cjoinby_fault) fault(post) testname("[rollback] string key, rows dropped")
_ctools_rb_case, cmd(cjoinby id skey using `rb_using', update) ///
    global(CTOOLS_cjoinby_fault) fault(store) rc(920) testname("[rollback] two keys, write fails")
* a join that moves no master observation still changes storage types,
* labels and (with update) the shared variables
quietly {
    use `rb_using', clear
    bysort id (uid): keep if _n == 1
    tempfile rb_one
    save `rb_one'
    use `rb_copy', clear
}
foreach fault in schema store post strl {
    local rc = cond("`fault'" == "store", 920, 1)
    _ctools_rb_case, cmd(cjoinby id using `rb_one', unmatched(both) update replace) ///
        global(CTOOLS_cjoinby_fault) fault(`fault') rc(`rc') ///
        testname("[rollback] m:1 join, `fault'")
}
* a failure after the plan and before the first change (hook record)
global CTOOLS_cjoinby_fault record
capture cjoinby id using `rb_using', unmatched(both) update
local rc = _rc
global CTOOLS_cjoinby_fault
_ctools_rb_state
local after `"`r(state)'"'
quietly use `rb_copy', clear
_ctools_rb_state
local before `"`r(state)'"'
if `rc' == 1 & `"`before'"' == `"`after'"' test_pass "[rollback] failure before the first change"
else test_fail "[rollback] failure before the first change" "r(`rc') or the data changed"
* an error before any change: nothing to undo, nothing left behind
capture cjoinby id using `rb_using', _merge(uid)
local rc = _rc
capture cf _all using `rb_copy', all
local left : all globals "CTOOLS_cjoinby*"
if `rc' == 110 & !_rc & "`left'" == "" test_pass "[rollback] error before any change"
else test_fail "[rollback] error before any change" "r(`rc'), cf r(`=_rc'), globals `left'"

* ------------------------------------------------------------------------
print_section "Without frames (Stata 14 and 15)"
* ------------------------------------------------------------------------
* Before Stata 16 the using data can only be read in place of the preserved
* master. The hook $CTOOLS_cjoinby_noframes runs that path here.
global CTOOLS_cjoinby_noframes 1
foreach keytype in int string two {
    local keys id
    if "`keytype'" == "two" local keys g id
    _cjoinby_data `m' 60 1101 m `keytype'
    _cjoinby_data `u' 45 2202 u `keytype'
    quietly {
        use `m', clear
        sort `keys' mid
        save `ref', replace
    }
    foreach spec in "" "unmatched(both) update replace" "unmatched(using) nolabel" {
        _cjoinby_parity, master(`m') using(`u') keys(`keys') opts(`spec') ///
            testname("no frames, `keytype' keys, `spec'")
        _cjoinby_parity, master(`ref') using(`u') keys(`keys') opts(`spec') ///
            testname("no frames, sorted master, `keytype' keys, `spec'")
    }
}
_cjoinby_parity, master(`m') using(`u') keys(mx) testname("no frames, key not in the using data (r(111))")
quietly use `rb_copy', clear
foreach fault in schema store post strl {
    local rc = cond("`fault'" == "store", 920, 1)
    _ctools_rb_case, cmd(cjoinby id using `rb_using', unmatched(both) update replace) ///
        global(CTOOLS_cjoinby_fault) fault(`fault') rc(`rc') modes(default) ///
        testname("[rollback] no frames, `fault'")
}
global CTOOLS_cjoinby_noframes

* ------------------------------------------------------------------------
print_section "e(sample)"
* ------------------------------------------------------------------------
* joinby moves e(sample) with each master observation (0 on using-only
* ones); stored samples (_est_*) are joined like any master variable.
do "validation/support/esample_helpers.do"
clear
quietly set obs 9
gen int k = mod(_n, 4) + 1
gen double v = 10 * _n
tempfile es_using
quietly save `es_using'
global CTOOLS_esample_using `es_using'
capture program drop _cjoinby_esample_setup
program define _cjoinby_esample_setup
    clear
    quietly set obs 12
    gen long id = _n
    gen int k = mod(_n, 5)
    gen double x = mod(_n * 7, 12)
    gen double y = x + mod(_n * 5, 7)
    gen byte g = mod(_n, 3) != 0
    quietly regress y x if !g
    estimates store m2
    quietly regress y x if g
    estimates store m1
end
* joinby leaves e(sample) with one pairing of each master observation
* (expand's copies are new observations), and which one depends on the sort
* seed; cjoinby leaves it with the first. Compared per master observation;
* the samples of stored estimates are compared row by row.
capture program drop _cjoinby_esample_side
program define _cjoinby_esample_side
    quietly gen byte es_after = e(sample)
    quietly predict double es_p if e(sample)
    quietly estimates restore m2
    quietly gen byte es_m2 = e(sample)
    quietly estimates restore m1
    quietly gen byte es_restored = e(sample)
    sort id v
end
tempfile es_ref es_ref_sum
quietly {
    _cjoinby_esample_setup
    joinby k using $CTOOLS_esample_using, unmatched(both)
    local st_native "functions=[`: e(functions)'] name=[`e(_estimates_name)'] cmd=[`e(cmd)'] changed=`c(changed)'"
    _cjoinby_esample_side
    preserve
    collapse (sum) es_after (count) es_p, by(id)
    sort id
    save `es_ref_sum'
    restore
    drop es_after es_p
    save `es_ref'
}
foreach mode in "" "snapshot" "fast" {
    quietly {
        _cjoinby_esample_setup
        cjoinby k using $CTOOLS_esample_using, unmatched(both) `mode'
        local st_ctools "functions=[`: e(functions)'] name=[`e(_estimates_name)'] cmd=[`e(cmd)'] changed=`c(changed)'"
        _cjoinby_esample_side
        capture assert es_after == 0 if id == id[_n - 1] & _n > 1
        local rc_first = _rc
        preserve
        collapse (sum) es_after (count) es_p, by(id)
        sort id
        capture cf _all using `es_ref_sum', all
        local rc_sums = _rc
        restore
        drop es_after es_p
        capture cf _all using `es_ref', all
        local rc_rows = _rc
    }
    if `"`st_native'"' == `"`st_ctools'"' & !`rc_first' & !`rc_sums' & !`rc_rows' {
        test_pass "e(sample) after cjoinby (`mode')"
    }
    else {
        test_fail "e(sample) after cjoinby (`mode')" ///
            `"native `st_native' vs ctools `st_ctools'; first r(`rc_first'), sums r(`rc_sums'), rows r(`rc_rows')"'
    }
}
foreach fault in schema post {
    global CTOOLS_cjoinby_fault `fault'
    _ctools_esample_case, setup(_cjoinby_esample_setup) native(error 1) ///
        ctools(cjoinby k using $CTOOLS_esample_using, unmatched(both)) ///
        testname("e(sample) after a failed cjoinby (`fault')")
    global CTOOLS_cjoinby_fault
}
global CTOOLS_esample_using

print_summary "cjoinby"
global CTOOLS_COMPONENT_COMPLETE "cjoinby"
    global CTOOLS_CASE_COMPLETE "cjoinby"
    exit
}

if "`__case'" == "cjoinby_extra" {
* Extra reference tests against Stata's joinby: key edge values, update and
* replace cell by cell, names, paths, frames, settings and error states.
do "validation/support/validate_setup.do"
capture confirm file "validation/support/rollback_helpers.do"
if !_rc do "validation/support/rollback_helpers.do"
else do "rollback_helpers.do"

* _cjx_parity: as _cjoinby_parity (case cjoinby), and when both commands fail
* also the master must be as it was (state, observations, no leftover
* globals or frames); rc() checks the expected return code.
capture program drop _cjx_parity
program define _cjx_parity
    syntax, master(string) using(string) testname(string) [keys(string) ///
        opts(string asis) copts(string asis) rc(integer -1)]
    tempfile ref
    quietly use `"`master'"', clear
    capture joinby `keys' using `"`using'"', `opts'
    local rc_native = _rc
    if !`rc_native' {
        _ctools_rb_state
        local st_native `"`r(state)'"'
        local all
        foreach v of varlist _all {
            local all `all' `v'
        }
        sort `all'
        quietly save `"`ref'"'
    }
    quietly use `"`master'"', clear
    _ctools_rb_state
    local st_before `"`r(state)'"'
    quietly frames dir
    local frames_before `"`r(frames)'"'
    capture cjoinby `keys' using `"`using'"', `opts' `copts'
    local rc_ctools = _rc
    quietly frames dir
    local frames_after `"`r(frames)'"'
    local left : all globals "CTOOLS_cjoinby*"
    if `rc_native' != `rc_ctools' {
        test_fail "`testname'" "joinby r(`rc_native'), cjoinby r(`rc_ctools')"
        exit
    }
    if `rc' >= 0 & `rc_native' != `rc' {
        test_fail "`testname'" "both r(`rc_native'), expected r(`rc')"
        exit
    }
    if `"`frames_before'"' != `"`frames_after'"' | "`left'" != "" {
        test_fail "`testname'" "left behind: frames [`frames_after'] globals [`left']"
        exit
    }
    if `rc_native' {
        _ctools_rb_state
        local st_after `"`r(state)'"'
        mata: st_local("same", strofreal(st_local("st_before") == st_local("st_after")))
        capture cf _all using `"`master'"', all
        if !`same' | _rc {
            test_fail "`testname'" "the master changed after r(`rc_ctools')"
        }
        else test_pass "`testname'"
        exit
    }
    _ctools_rb_state
    local st_ctools `"`r(state)'"'
    mata: st_local("same", strofreal(st_local("st_native") == st_local("st_ctools")))
    if !`same' {
        mata: _cjx_first_difference(st_local("st_native"), st_local("st_ctools"))
        test_fail "`testname'" `"state differs: joinby [`macval(d_native)'] cjoinby [`macval(d_ctools)']"'
        exit
    }
    local sortedby : sortedby
    if "`sortedby'" != "" & _N > 0 {
        tempvar pos
        quietly generate long `pos' = _n
        sort `sortedby', stable
        capture assert `pos' == _n
        if _rc {
            test_fail "`testname'" "the data are not in the order of the sort marker (`sortedby')"
            exit
        }
        drop `pos'
    }
    sort `all'
    capture cf _all using `"`ref'"', all
    if _rc test_fail "`testname'" "observations differ from joinby (cf r(`=_rc'))"
    else test_pass "`testname'"
end

capture mata: mata drop _cjx_first_difference()
mata:
void _cjx_first_difference(string scalar a, string scalar b)
{
    real scalar i, n

    n = min((strlen(a), strlen(b)))
    for (i = 1; i <= n; i++) {
        if (substr(a, i, 1) != substr(b, i, 1)) break
    }
    st_local("d_native", substr(a, max((1, i - 60)), 140))
    st_local("d_ctools", substr(b, max((1, i - 60)), 140))
}
end

* _cjx_rows <list of variables>: the observations as "v1:v2:..." words in
* local rows of the caller, in the current order
capture program drop _cjx_rows
program define _cjx_rows
    syntax varlist
    tempvar row
    quietly gen strL `row' = ""
    foreach v of local varlist {
        capture confirm string variable `v'
        if _rc quietly replace `row' = `row' + ":" + string(`v')
        else quietly replace `row' = `row' + ":" + `v'
    }
    local rows
    forvalues i = 1/`=_N' {
        local rows `"`rows' `=substr(`row'[`i'], 2, .)'"'
    }
    c_local rows `"`=strtrim(`"`rows'"')'"'
end

tempfile m u ref out

* ------------------------------------------------------------------------
print_section "Extra: string keys"
* ------------------------------------------------------------------------
* Bytes decide: case, leading and trailing blanks, composed and decomposed
* accents, quotes, an invalid UTF-8 byte and the empty string are each a
* key of their own.
local e1 = uchar(233)
local e2 = "e" + uchar(769)
quietly {
    clear
    set obs 13
    gen long mid = _n
    gen str20 k = ""
    replace k = "a" in 1
    replace k = "A" in 2
    replace k = "a " in 3
    replace k = " a" in 4
    replace k = "" in 5
    replace k = "`e1'" in 6
    replace k = "`e2'" in 7
    replace k = `"q"uote"' in 8
    replace k = "it's" in 9
    replace k = uchar(26085) + uchar(26412) in 10
    replace k = "z" + char(200) in 11
    replace k = "a" in 12
    replace k = 18 * "w" in 13
    gen double mv = _n
    save `m', replace
    clear
    set obs 14
    gen long uid = _n
    gen str2045 k = ""
    replace k = "a" in 1
    replace k = "a" in 2
    replace k = "a  " in 3
    replace k = " a" in 4
    replace k = "" in 5
    replace k = "`e2'" in 6
    replace k = `"q"uote"' in 7
    replace k = "IT'S" in 8
    replace k = uchar(26085) + uchar(26412) in 9
    replace k = "z" + char(200) in 10
    replace k = "z" in 11
    replace k = 18 * "w" + " " in 12
    replace k = 18 * "w" in 13
    replace k = "`e1'" in 14
    gen double mv = -_n
    gen str3 uonly = "u" + string(_n)
    save `u', replace
}
foreach unmatched in none both master using {
    foreach how in "" "update replace" {
        _cjx_parity, master(`m') using(`u') keys(k) opts(unmatched(`unmatched') `how') ///
            testname("string keys by bytes, unmatched(`unmatched') `how'")
    }
}
quietly {
    use `m', clear
    sort k
    save `ref', replace
}
_cjx_parity, master(`ref') using(`u') keys(k) opts(unmatched(both)) ///
    testname("string keys by bytes, master sorted by the key")
* the counts behind it: "a" twice in each file, eight other keys once in each
quietly {
    use `m', clear
    cjoinby k using `u', unmatched(both)
}
local ok = r(N_3) == 4 + 8 & r(N_1) == 3 & r(N_2) == 4
if `ok' test_pass "string keys by bytes, r(N_1) r(N_2) r(N_3)"
else test_fail "string keys by bytes, r(N_1) r(N_2) r(N_3)" "r(N_1)=`r(N_1)' r(N_2)=`r(N_2)' r(N_3)=`r(N_3)'"

* ------------------------------------------------------------------------
print_section "Extra: numeric keys of different storage types"
* ------------------------------------------------------------------------
* byte master key against a double using key with near and non-integer values
quietly {
    clear
    set obs 8
    gen byte id = .
    replace id = 1 in 1
    replace id = 2 in 2
    replace id = 100 in 3
    replace id = -127 in 4
    replace id = 0 in 5
    replace id = .a in 6
    replace id = .z in 7
    replace id = . in 8
    gen long mid = _n
    save `m', replace
    clear
    set obs 10
    gen double id = .
    replace id = 1 in 1
    replace id = 1 + 1e-10 in 2
    replace id = 2 - epsdouble() in 3
    replace id = 100 in 4
    replace id = -127 in 5
    replace id = -0 in 6
    replace id = .a in 7
    replace id = .z in 8
    replace id = 2.5 in 9
    replace id = . in 10
    gen long uid = _n
    save `u', replace
}
foreach how in "" "update" {
    _cjx_parity, master(`m') using(`u') keys(id) opts(unmatched(both) `how') ///
        testname("byte master key, double using key `how'")
}
_cjx_parity, master(`u') using(`m') keys(id) opts(unmatched(both)) ///
    testname("double master key, byte using key")
* A non-integer master key one rounding step from an integer using key must
* not match it, also when the using keys reach far below zero (the distance
* to the smallest using key then rounds to a whole number)
foreach spec in "byte -127 2-epsdouble()" "byte -127 1e-15" "int -32767 5+1e-12" ///
    "long -1000000 5+1e-11" "byte -1 2-epsdouble()" {
    gettoken ut spec : spec
    gettoken lo mval : spec
    quietly {
        clear
        set obs 4
        gen `ut' id = `lo'
        replace id = 0 in 2
        replace id = 2 in 3
        replace id = 5 in 4
        gen long uid = _n
        save `out', replace
        clear
        set obs 3
        gen double id = `mval'
        replace id = 5 in 2
        replace id = 2 - 2 * epsdouble() in 3
        gen long mid = _n
        save `ref', replace
    }
    _cjx_parity, master(`ref') using(`out') keys(id) opts(unmatched(both)) ///
        testname("near-integer double master key, `ut' using keys from `lo', `mval'")
}
* float precision: float(0.1) is not 0.1
quietly {
    clear
    set obs 4
    gen float id = .
    replace id = 0.1 in 1
    replace id = 16777216 in 2
    replace id = 1/3 in 3
    replace id = -0.1 in 4
    gen long mid = _n
    save `m', replace
    clear
    set obs 8
    gen double id = .
    replace id = 0.1 in 1
    replace id = float(0.1) in 2
    replace id = 16777217 in 3
    replace id = 16777216 in 4
    replace id = 1/3 in 5
    replace id = float(1/3) in 6
    replace id = float(-0.1) in 7
    replace id = -0.1 in 8
    gen long uid = _n
    save `u', replace
}
_cjx_parity, master(`m') using(`u') keys(id) opts(unmatched(both)) ///
    testname("float master key against float-rounded and exact doubles")
quietly {
    use `m', clear
    cjoinby id using `u'
    sort mid uid
    _cjx_rows mid uid
}
if `"`rows'"' == "1:2 2:4 3:6 4:7" test_pass "float keys match only float-rounded doubles"
else test_fail "float keys match only float-rounded doubles" `"pairs `rows'"'
* integer and long edge values, with and without an int key widening
quietly {
    clear
    set obs 6
    gen long id = .
    replace id = 2147483620 in 1
    replace id = -2147483647 in 2
    replace id = 0 in 3
    replace id = 32740 in 4
    replace id = -32767 in 5
    replace id = .b in 6
    gen long mid = _n
    save `m', replace
    clear
    set obs 6
    gen int id = .
    replace id = 32740 in 1
    replace id = -32767 in 2
    replace id = 0 in 3
    replace id = .b in 4
    replace id = 32740 in 5
    replace id = 1 in 6
    gen long uid = _n
    save `u', replace
}
_cjx_parity, master(`m') using(`u') keys(id) opts(unmatched(both)) ///
    testname("long master key edges against int using key")
_cjx_parity, master(`u') using(`m') keys(id) opts(unmatched(both) update) ///
    testname("int master key against long using key edges")

* Doubles beyond the integer range of a 64-bit table position: large and
* tiny magnitudes, values a few units of 2^53 apart, values around 2^63
quietly {
    clear
    set obs 18
    gen double id = .
    replace id = 1e18 in 1
    replace id = 1e18 + 128 in 2
    replace id = 1e18 + 256 in 3
    replace id = 2^53 in 4
    replace id = 2^53 + 2 in 5
    replace id = 2^63 in 6
    replace id = 2^63 - 1024 in 7
    replace id = -2^63 in 8
    replace id = -2^63 - 2048 in 9
    replace id = 1e19 in 10
    replace id = c(maxdouble) in 11
    replace id = c(mindouble) in 12
    replace id = c(smallestdouble) in 13
    replace id = -c(smallestdouble) in 14
    replace id = 1e-300 in 15
    replace id = 2^64 in 16
    replace id = 2^64 + 4096 in 17
    replace id = 0 in 18
    gen long mid = _n
    save `m', replace
    keep in 1/10
    replace mid = -mid
    rename mid uid
    expand 2 in 2
    replace id = 1e18 + 384 in 5
    replace id = 2^53 + 4 in 4
    set obs 14
    replace id = 2^64 in 12
    replace id = -0 in 13
    replace id = c(maxdouble) in 14
    replace uid = 100 + _n if missing(uid)
    save `u', replace
}
foreach unmatched in none both {
    _cjx_parity, master(`m') using(`u') keys(id) opts(unmatched(`unmatched')) ///
        testname("huge and tiny double keys, unmatched(`unmatched')")
}
* the same values in a narrow range only: a direct table of large integers
quietly {
    clear
    set obs 40
    gen double id = 4e18 + 1024 * mod(_n, 9)
    gen long mid = _n
    save `m', replace
    clear
    set obs 30
    gen double id = 4e18 + 1024 * mod(_n, 11)
    gen long uid = _n
    save `u', replace
}
_cjx_parity, master(`m') using(`u') keys(id) opts(unmatched(both)) ///
    testname("large integers 1024 apart")
quietly {
    clear
    set obs 40
    gen double id = 2^53 + mod(_n, 9)
    gen long mid = _n
    save `m', replace
    clear
    set obs 30
    gen double id = 2^53 + 2 * mod(_n, 11)
    gen long uid = _n
    save `u', replace
}
_cjx_parity, master(`m') using(`u') keys(id) opts(unmatched(both)) ///
    testname("integers past 2^53 (rounded on input)")

* Three long keys spanning their whole range: the combined range of the
* key tuple does not fit 64 bits
quietly {
    clear
    set seed 17001
    set obs 300
    foreach v in a b c {
        gen long `v' = cond(runiform() < .5, -2147483647, 2147483620)
        replace `v' = 0 if runiform() < .3
        replace `v' = . if runiform() < .05
    }
    gen long mid = _n
    save `m', replace
    clear
    set obs 200
    foreach v in a b c {
        gen long `v' = cond(runiform() < .5, -2147483647, 2147483620)
        replace `v' = 0 if runiform() < .3
        replace `v' = . if runiform() < .05
    }
    gen long uid = _n
    save `u', replace
}
_cjx_parity, master(`m') using(`u') keys(a b c) opts(unmatched(both)) ///
    testname("three long keys over their whole range")
* the same in three byte and int keys whose tuple range is small, and a
* tuple whose first key covers the range and whose last is a double
quietly {
    use `m', clear
    replace a = mod(abs(a), 7)
    recast byte a, force
    save `ref', replace
    use `u', clear
    replace a = mod(abs(a), 7)
    recast double a
    replace c = c + .5 if mod(uid, 3) == 0
    save `out', replace
}
_cjx_parity, master(`ref') using(`out') keys(a b c) opts(unmatched(both)) ///
    testname("byte, long and long keys against double, long and double")


* Every missing value in byte, float and double keys
quietly {
    clear
    set obs 27
    gen byte id = .
    forvalues j = 2/27 {
        replace id = real("." + char(95 + `j')) in `j'
    }
    gen long mid = _n
    expand 2 in 1/5
    set obs `=_N + 1'
    replace id = 7 in l
    save `m', replace
    clear
    set obs 27
    gen double id = .
    forvalues j = 2/27 {
        replace id = real("." + char(123 - `j' + 1)) in `j'
    }
    gen long uid = _n
    expand 3 if id == .z | id == .
    set obs `=_N + 1'
    replace id = 7.5 in l
    replace uid = 99 in l
    save `u', replace
}
foreach how in "" "update" {
    _cjx_parity, master(`m') using(`u') keys(id) opts(unmatched(both) `how') ///
        testname("every missing value, byte master key, double using key `how'")
}
quietly {
    use `m', clear
    recast float id
    save `ref', replace
}
_cjx_parity, master(`ref') using(`u') keys(id) opts(unmatched(both)) ///
    testname("every missing value, float master key")
quietly {
    use `ref', clear
    replace id = id + .25 if !missing(id)
    save `ref', replace
}
_cjx_parity, master(`ref') using(`u') keys(id) opts(unmatched(both)) ///
    testname("every missing value with non-integer keys")
quietly {
    use `m', clear
    cjoinby id using `u', unmatched(both)
}
* . 2 x 3, .z 1 x 3, .a-.d 2 x 1, .e-.y 1 x 1; 7 and 7.5 unmatched
local ok = r(N_3) == 6 + 3 + 8 + 21 & r(N_1) == 1 & r(N_2) == 1
if `ok' test_pass "every missing value, pairing counts"
else test_fail "every missing value, pairing counts" "r(N_1)=`r(N_1)' r(N_2)=`r(N_2)' r(N_3)=`r(N_3)'"

* ------------------------------------------------------------------------
print_section "Extra: sizes of many-to-many joins"
* ------------------------------------------------------------------------
quietly {
    clear
    set obs 13
    gen int k = cond(_n <= 7, 1, cond(_n <= 10, 2, cond(_n <= 12, ., 9)))
    gen long mid = _n
    save `m', replace
    clear
    set obs 17
    gen int k = cond(_n <= 5, 1, cond(_n <= 9, 2, cond(_n <= 15, ., 3)))
    gen long uid = _n
    save `u', replace
}
local okall 1
foreach spec in "none 59 0 0" "both 62 1 2" "master 60 1 0" "using 61 0 2" {
    gettoken unm spec : spec
    gettoken n spec : spec
    gettoken n1 spec : spec
    gettoken n2 spec : spec
    quietly use `m', clear
    quietly cjoinby k using `u', unmatched(`unm')
    local ok = _N == `n' & r(N) == `n' & r(N_1) == `n1' & r(N_2) == `n2' & r(N_3) == 59
    if `ok' {
        quietly count if k == 1
        local ok = r(N) == 35
        quietly count if missing(k)
        local ok = `ok' & r(N) == 12
    }
    if `ok' test_pass "group sizes multiply, unmatched(`unm')"
    else test_fail "group sizes multiply, unmatched(`unm')" "_N=`=_N', expected `n'"
}
quietly {
    clear
    set seed 17002
    set obs 401
    gen int k = cond(_n <= 400, 1, 2)
    gen long mid = _n
    gen double mx = runiform()
    save `m', replace
    clear
    set obs 302
    gen int k = cond(_n <= 300, 1, cond(_n == 301, 2, 5))
    gen long uid = _n
    gen str4 us = "u" + string(mod(_n, 13))
    save `u', replace
}
_cjx_parity, master(`m') using(`u') keys(k) opts(unmatched(both)) ///
    testname("400 x 300 pairings of one key")
quietly {
    use `m', clear
    cjoinby k using `u'
}
if _N == 120001 & r(N_3) == 120001 test_pass "400 x 300 + 1 x 1 = 120,001 observations"
else test_fail "400 x 300 + 1 x 1 = 120,001 observations" "_N=`=_N'"

* ------------------------------------------------------------------------
print_section "Extra: update and replace, cell by cell"
* ------------------------------------------------------------------------
quietly {
    clear
    input long mid int id double v str2 s
    1 1 .  ""
    2 1 .a "m"
    3 2 5  ""
    4 3 6  "m"
    5 4 .  ""
    6 5 7  "m"
    7 6 .  ""
    8 8 .a ""
    9 9 .  "m"
    end
    save `m', replace
    clear
    input int id double v str2 s
    1 10 "u"
    2 .  "u"
    3 .b ""
    4 40 ""
    5 8  "x"
    7 70 "y"
    8 .b "z"
    9 .c ""
    end
    save `u', replace
}
foreach how in "" "update" "update replace" {
    _cjx_parity, master(`m') using(`u') keys(id) opts(unmatched(both) `how') ///
        testname("update cells, parity `how'")
}
local exp_none "1:1:.: 2:1:.a:m 3:2:5: 4:3:6:m 5:4:.: 6:5:7:m 7:6:.: .:7:.:"
local exp_update "1:1:10:u 2:1:10:m 3:2:5:u 4:3:6:m 5:4:40: 6:5:7:m 7:6:.: .:7:70:y"
local exp_replace "1:1:10:u 2:1:10:u 3:2:5:u 4:3:6:m 5:4:40: 6:5:8:x 7:6:.: .:7:70:y"
foreach how in none update replace {
    local o
    if "`how'" == "update" local o update
    if "`how'" == "replace" local o update replace
    quietly {
        use `m', clear
        cjoinby id using `u', unmatched(both) `o'
        drop if id >= 8
        sort mid
        _cjx_rows mid id v s
    }
    if `"`rows'"' == `"`exp_`how''"' test_pass "update cells, expected values (`how')"
    else test_fail "update cells, expected values (`how')" `"got `rows'"'
}
* _merge stays 1/2/3 and __MERGE has the update wording
quietly {
    use `m', clear
    cjoinby id using `u', unmatched(both) update replace
    count if !inlist(_merge, 1, 2, 3)
    local bad = r(N)
    local l4 : label __MERGE 4
    local lv : value label _merge
}
if `bad' == 0 & "`l4'" == "in both, missing in master data updated" & "`lv'" == "__MERGE" {
    test_pass "update: _merge codes and __MERGE"
}
else test_fail "update: _merge codes and __MERGE" "bad=`bad' label 4 [`l4'] value label [`lv']"

* ------------------------------------------------------------------------
print_section "Extra: option parsing"
* ------------------------------------------------------------------------
quietly {
    clear
    set seed 17004
    set obs 30
    gen long mid = _n
    gen int id = floor(runiform() * 8)
    gen double sh = cond(runiform() < .4, ., _n)
    save `m', replace
    clear
    set obs 25
    gen long uid = _n
    gen int id = floor(runiform() * 10)
    gen double sh = -_n
    gen byte _merge = 9
    save `ref', replace
    drop _merge
    save `u', replace
}
foreach o in "unmatched(BOTH)" "unmatched(Ma)" "unmatched()" "unmatched(us)" ///
    "unmatched(n)" "unm(u) nol" "unmatched(bothx)" "unmatched(both master)" ///
    "upd" "update rep" "replace update" "nolab update" ///
    "_merge(1abc) unmatched(both)" "_merge(a b)" "_merge(id) unmatched(both)" ///
    "_merge(src) unmatched(master)" "_merge(_merge)" {
    _cjx_parity, master(`m') using(`u') keys(id) opts(`o') testname("options: `o'")
}
* a _merge variable in the using data
_cjx_parity, master(`m') using(`ref') keys(id) opts(unmatched(both)) rc(110) ///
    testname("_merge in the using data, unmatched(both) (r(110))")
_cjx_parity, master(`m') using(`ref') keys(id) ///
    testname("_merge in the using data without unmatched() is a using variable")
_cjx_parity, master(`m') using(`ref') keys(id) opts(_merge(src) unmatched(both)) ///
    testname("_merge in the using data, _merge(src)")

* ------------------------------------------------------------------------
print_section "Extra: variable names and varlists"
* ------------------------------------------------------------------------
quietly {
    clear
    set seed 17005
    set obs 40
    gen long mid = _n
    gen int kfirm = floor(runiform() * 6)
    gen byte kyear = floor(runiform() * 3)
    gen double sh = runiform()
    save `m', replace
    clear
    set obs 30
    gen long uid = _n
    gen int kfirm = floor(runiform() * 7)
    gen byte kyear = floor(runiform() * 3)
    gen double sh = -_n
    save `u', replace
}
_cjx_parity, master(`m') using(`u') keys(kfi kye) opts(unmatched(both)) ///
    testname("abbreviated key names")
_cjx_parity, master(`m') using(`u') keys(k*) opts(unmatched(both)) ///
    testname("wildcard key varlist")
_cjx_parity, master(`m') using(`u') keys(kfirm-kyear) opts(unmatched(both)) ///
    testname("key varlist range")
_cjx_parity, master(`m') using(`u') keys(kyear kfirm) opts(unmatched(both) update) ///
    testname("keys in the other order")
_cjx_parity, master(`m') using(`u') keys(kfirm kfirm) opts(unmatched(both)) ///
    testname("a key named twice")
set varabbrev off
_cjx_parity, master(`m') using(`u') keys(kfi) testname("varabbrev off, abbreviated key (r(111))")
_cjx_parity, master(`m') using(`u') keys(kfirm kyear) opts(unmatched(both) update replace) ///
    testname("varabbrev off, full names")
_cjx_parity, master(`m') using(`u') opts(unmatched(using)) ///
    testname("varabbrev off, common variables")
set varabbrev on

* Unicode and 32-character names
local vk = "cl" + uchar(233) + "_" + uchar(26085)
local long32 = "a" + 31 * "b"
local long32b = "a" + 30 * "b" + "c"
quietly {
    use `m', clear
    rename kfirm `vk'
    rename sh `long32'
    gen float `long32b' = mid
    save `ref', replace
    use `u', clear
    rename kfirm `vk'
    rename sh `long32'
    gen str3 `long32b' = "u"
    gen int z`vk' = uid
    save `out', replace
}
_cjx_parity, master(`ref') using(`out') keys(`vk' kyear) opts(unmatched(both) update) ///
    testname("Unicode and 32-character variable names")

* Variables of the using (or master) data named like the next temporary
* names: cjoinby's own temporary variables must not collide with them.
capture program drop _cjx_nexttemp
program define _cjx_nexttemp
    tempvar t
    c_local next `t'
end
capture mata: mata drop _cjx_b36()
mata:
string scalar _cjx_b36(real scalar n0)
{
    string scalar d, s
    real scalar i, n

    d = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    s = ""
    n = n0
    for (i = 1; i <= 6; i++) {
        s = substr(d, mod(n, 36) + 1, 1) + s
        n = floor(n / 36)
    }
    return("__" + s)
}
real scalar _cjx_b36val(string scalar t)
{
    string scalar d
    real scalar i, n

    d = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    n = 0
    for (i = 3; i <= strlen(t); i++) n = n * 36 + strpos(d, substr(t, i, 1)) - 1
    return(n)
}
end
quietly {
    clear
    set seed 17006
    set obs 50
    gen long mid = _n
    gen int id = floor(runiform() * 10)
    gen double sh = cond(runiform() < .5, ., _n)
    gen strL ms = "m" + string(_n)
    save `m', replace
    clear
    set obs 40
    gen long uid = _n
    gen int id = floor(runiform() * 12)
    gen double sh = -_n
    gen strL us = "u" + string(_n)
    save `u', replace
    use `m', clear
    cjoinby id using `u', unmatched(both) update
    sort mid uid
    save `ref', replace
}
foreach side in using master {
    quietly {
        use `m', clear
        _cjx_nexttemp
        mata: st_local("n0", strofreal(_cjx_b36val("`next'")))
        if "`side'" == "using" use `u', clear
        forvalues j = 1/200 {
            mata: st_local("nm", _cjx_b36(`n0' + `j'))
            gen byte `nm' = 7
        }
        if "`side'" == "using" {
            save `out', replace
            use `m', clear
        }
        else save `out', replace
        * (cf would collide with the names itself: compare states instead)
        _ctools_rb_state
        local st_before `"`r(state)'"'
        if "`side'" == "using" {
            capture cjoinby id using `out', unmatched(both) update
        }
        else {
            capture cjoinby id using `u', unmatched(both) update
        }
        local rc = _rc
    }
    * joinby fails here too (its own temporary names collide), so a failure
    * is accepted when it leaves the master as it was
    if `rc' {
        _ctools_rb_state
        local st_after `"`r(state)'"'
        mata: st_local("same", strofreal(st_local("st_before") == st_local("st_after")))
        if `rc' == 110 & `same' test_pass "`side' variables named like temporary variables (r(110), master unchanged)"
        else test_fail "`side' variables named like temporary variables" "r(`rc'), master unchanged `same'"
        continue
    }
    local ok = !`rc'
    if `ok' {
        mata: st_local("first", _cjx_b36(`n0' + 1))
        mata: st_local("last", _cjx_b36(`n0' + 200))
        local mark = cond("`side'" == "using", 1, 2)
        capture assert `first' == cond(_merge == `mark', ., 7) & `last' == `first'
        local ok = !_rc
        quietly drop `first'-`last'
        sort mid uid
        capture cf _all using `ref', all
        local ok = `ok' & !_rc
    }
    if `ok' test_pass "`side' variables named like temporary variables"
    else test_fail "`side' variables named like temporary variables" "r(`rc'), the names from `first' on"
}

* ------------------------------------------------------------------------
print_section "Extra: labels, notes and characteristics of both files"
* ------------------------------------------------------------------------
quietly {
    clear
    set seed 17007
    set obs 30
    gen long mid = _n
    gen int id = floor(runiform() * 6)
    gen double sh = cond(runiform() < .5, ., floor(runiform() * 3))
    gen double mx = runiform()
    note id: the same text
    note id: master only note on id
    note: the same dataset note
    char id[mchar] master id characteristic
    label define lsh 0 "master zero"
    label values mx lsh
    save `m', replace
    clear
    set obs 25
    gen long uid = _n
    gen int id = floor(runiform() * 7)
    gen double sh = floor(runiform() * 3)
    gen str5 ux = "u" + string(_n)
    label define lsh 0 "using zero" 1 "using one"
    label define lux 1 `"with "quotes" and `=uchar(233)'"'
    label values sh lsh
    label values uid lux
    label variable sh `"using "quoted" label `=uchar(26085)'"'
    label variable id "using key label"
    forvalues j = 1/12 {
        note id: using note `j' on id
    }
    note id: the same text
    forvalues j = 1/11 {
        note: using dataset note `j'
    }
    note: the same dataset note
    char ux[long] `=3000 * "c"'
    char _dta[ulong] `=2500 * "d"'
    char id[mchar] using id characteristic
    char id[uchar] using-only id characteristic
    label data "using label"
    save `u', replace
}
foreach o in "unmatched(both)" "unmatched(both) update" "unmatched(both) nolabel update" "nolabel" {
    _cjx_parity, master(`m') using(`u') keys(id) opts(`o') ///
        testname("labels, notes, characteristics, `o'")
}
quietly {
    use `m', clear
    label data ""
    notes drop _all
    save `ref', replace
}
_cjx_parity, master(`ref') using(`u') keys(id) opts(unmatched(both) update) ///
    testname("labels, notes, characteristics, master without notes or data label")

* ------------------------------------------------------------------------
print_section "Extra: file names"
* ------------------------------------------------------------------------
quietly {
    clear
    set seed 17008
    set obs 20
    gen long mid = _n
    gen int id = floor(runiform() * 5)
    save `m', replace
    clear
    set obs 15
    gen long uid = _n
    gen int id = floor(runiform() * 6)
    save `u', replace
}
local dotdir "`c(tmpdir)'/cjx dir.v1 extra"
capture mkdir "`dotdir'"
quietly {
    use `u', clear
    save "`dotdir'/using data.dta", replace
    save "`dotdir'/it's using.dta", replace
    save "`dotdir'/plain", replace
}
foreach f in "`dotdir'/using data" "`dotdir'/using data.dta" "`dotdir'/it's using" ///
    "`dotdir'/plain" "`dotdir'/absent file" {
    local short = substr(`"`f'"', strlen("`dotdir'") + 2, .)
    _cjx_parity, master(`m') using(`"`f'"') keys(id) opts(unmatched(both)) ///
        testname("using file `short' in a directory with dots and spaces")
}
local back = subinstr("`dotdir'/using data", "/", "\", .)
_cjx_parity, master(`m') using(`"`back'"') keys(id) opts(unmatched(both)) ///
    testname("using file with backslashes")
* a file that is not a Stata dataset
tempname fh
quietly file open `fh' using "`dotdir'/not data.dta", write text replace
file write `fh' "id,uid" _n "1,2" _n
file close `fh'
_cjx_parity, master(`m') using(`"`dotdir'/not data.dta"') keys(id) ///
    testname("using file that is not a dataset")
capture erase "`dotdir'/using data.dta"
capture erase "`dotdir'/it's using.dta"
capture erase "`dotdir'/plain.dta"
capture erase "`dotdir'/not data.dta"
capture rmdir "`dotdir'"

* ------------------------------------------------------------------------
print_section "Extra: frames"
* ------------------------------------------------------------------------
* Joined in a frame other than default; the default frame stays as it was
quietly {
    clear
    set seed 17009
    set obs 30
    gen long mid = _n
    gen int id = floor(runiform() * 6)
    gen strL ms = "m" + string(_n)
    save `m', replace
    clear
    set obs 25
    gen long uid = _n
    gen int id = floor(runiform() * 7)
    gen strL us = "u" + string(_n)
    save `u', replace
    use `m', clear
    joinby id using `u', unmatched(both)
    sort mid uid
    save `ref', replace
    sysuse auto, clear
    datasignature
    local sig_default `r(datasignature)'
}
capture frame drop cjx_work
frame create cjx_work
frame change cjx_work
quietly use `m', clear
capture cjoinby id using `u', unmatched(both)
local rc = _rc
local fr `c(frame)'
sort mid uid
capture cf _all using `ref', all
local rc_cf = _rc
frame change default
quietly datasignature
if !`rc' & !`rc_cf' & "`fr'" == "cjx_work" & "`r(datasignature)'" == "`sig_default'" & c(frame) == "default" {
    test_pass "join in a frame other than default"
}
else test_fail "join in a frame other than default" "r(`rc'), cf r(`rc_cf'), frame `fr'"
frame change cjx_work
quietly use `m', clear
capture cjoinby id using `u', _merge(uid)
local rc = _rc
capture cf _all using `m', all
local rc_cf = _rc
frame change default
quietly frames dir
local frs `r(frames)'
local frs : list sort frs
if `rc' == 110 & !`rc_cf' & "`frs'" == "cjx_work default" {
    test_pass "failed join in a frame other than default"
}
else test_fail "failed join in a frame other than default" "r(`rc'), cf r(`rc_cf'), frames `r(frames)'"
frame drop cjx_work

* ------------------------------------------------------------------------
print_section "Extra: errors leave the master unchanged"
* ------------------------------------------------------------------------
quietly {
    clear
    set seed 17010
    set obs 20
    gen long mid = _n
    gen int id = floor(runiform() * 5)
    gen str2 g = "g" + string(mod(_n, 3))
    gen double sh = _n
    sort id
    save `m', replace
    clear
    set obs 15
    gen long uid = _n
    gen int id = floor(runiform() * 6)
    gen str2 g = "g" + string(mod(_n, 2))
    gen double sh = -_n
    save `u', replace
    use `u', clear
    tostring id, replace
    save `ref', replace
    use `u', clear
    destring g, replace ignore("g")
    save `out', replace
}
_cjx_parity, master(`m') using(`u') keys(mid) rc(111) testname("error: key not in the using data")
_cjx_parity, master(`m') using(`ref') keys(id) rc(106) ///
    testname("error: numeric master key, string using key")
_cjx_parity, master(`m') using(`out') keys(id g) rc(106) ///
    testname("error: second key string in master, numeric in using")
_cjx_parity, master(`m') using(`ref') keys(g id) opts(unmatched(both) update replace) rc(106) ///
    testname("error: type mismatch with update replace")
_cjx_parity, master(`m') using(`u') keys(id nosuch) rc(111) testname("error: key not in the master")
* if/in are not allowed
quietly use `m', clear
capture joinby id using `u' if mid > 3
local rc_native = _rc
quietly use `m', clear
capture cjoinby id using `u' if mid > 3
local rc_ctools = _rc
capture cf _all using `m', all
if `rc_native' == `rc_ctools' & `rc_native' & !_rc test_pass "error: if qualifier (r(`rc_native'))"
else test_fail "error: if qualifier" "joinby r(`rc_native'), cjoinby r(`rc_ctools'), cf r(`=_rc')"
* no data in memory
clear
capture joinby id using `u'
local rc_native = _rc
clear
capture cjoinby id using `u'
local rc_ctools = _rc
if `rc_native' == `rc_ctools' & `rc_native' test_pass "error: no variables in memory, key named (r(`rc_native'))"
else test_fail "error: no variables in memory, key named" "joinby r(`rc_native'), cjoinby r(`rc_ctools')"
clear
capture joinby using `u'
local rc_native = _rc
clear
capture cjoinby using `u'
local rc_ctools = _rc
if `rc_native' == `rc_ctools' & `rc_native' test_pass "error: no variables in memory, common variables (r(`rc_native'))"
else test_fail "error: no variables in memory, common variables" "joinby r(`rc_native'), cjoinby r(`rc_ctools')"

* ------------------------------------------------------------------------
print_section "Extra: tiny data and self-joins"
* ------------------------------------------------------------------------
quietly {
    clear
    set obs 1
    gen int id = 3
    gen double mx = 1.5
    save `m', replace
    clear
    set obs 1
    gen double id = 3
    gen str1 us = "u"
    save `u', replace
    replace id = 4
    save `ref', replace
}
foreach unmatched in none both master using {
    _cjx_parity, master(`m') using(`u') keys(id) opts(unmatched(`unmatched')) ///
        testname("one observation each, matched, unmatched(`unmatched')")
    _cjx_parity, master(`m') using(`ref') keys(id) opts(unmatched(`unmatched')) ///
        testname("one observation each, unmatched keys, unmatched(`unmatched')")
}
quietly {
    use `m', clear
    keep id
    expand 3
    save `out', replace
    use `u', clear
    keep id
    expand 2
    save `ref', replace
}
_cjx_parity, master(`out') using(`ref') keys(id) opts(unmatched(both)) ///
    testname("only the key in each file")
_cjx_parity, master(`out') using(`ref') testname("only the key in each file, common variables")
quietly {
    clear
    set seed 17011
    set obs 40
    gen int id = floor(runiform() * 8)
    gen double x = cond(runiform() < .3, ., runiform())
    gen str3 s = cond(runiform() < .3, "", "s" + string(_n))
    gen strL l = "l" + string(mod(_n, 5))
    label variable x "x label"
    note x: self note
    save `m', replace
}
_cjx_parity, master(`m') using(`m') keys(id) opts(unmatched(both) update replace) ///
    testname("self-join on a key, update replace")
_cjx_parity, master(`m') using(`m') opts(unmatched(both)) ///
    testname("self-join on all variables, a strL among them")
quietly {
    use `m', clear
    drop l
    save `m', replace
}
_cjx_parity, master(`m') using(`m') opts(unmatched(both)) ///
    testname("self-join on all variables")

* ------------------------------------------------------------------------
print_section "Extra: master observations that stay in place"
* ------------------------------------------------------------------------
* The first observations match once (they do not move) while the key and
* shared variables widen; later ones expand.
quietly {
    clear
    set obs 300
    gen byte id = mod(_n - 1, 100)
    replace id = -5 in 241/250
    replace id = 100 in 251/300
    gen str3 s = "m" + string(mod(_n, 10))
    gen long sh = cond(mod(_n, 4) == 0, ., _n)
    gen long mid = _n
    save `m', replace
    clear
    set obs 104
    gen double id = _n - 1 + .5 * (_n > 101)
    replace id = 100 in 102/104
    gen str20 s = "using " + string(_n)
    gen float sh = _n + .25
    gen long uid = _n
    save `u', replace
}
foreach unmatched in master none both {
    foreach how in "" "update" "update replace" {
        _cjx_parity, master(`m') using(`u') keys(id) opts(unmatched(`unmatched') `how') ///
            testname("first observations stay, types widen, unmatched(`unmatched') `how'")
    }
}

* ------------------------------------------------------------------------
print_section "Extra: large joins"
* ------------------------------------------------------------------------
quietly {
    clear
    set seed 17012
    set obs 100000
    gen long uid = _n
    gen str14 id = "key" + string(floor(runiform() * 50000), "%07.0f")
    replace id = "" if runiform() < .01
    gen double sh = cond(runiform() < .3, ., runiform())
    gen str6 shs = cond(runiform() < .3, "", "u" + string(mod(_n, 991)))
    gen float ux = runiform()
    save `u', replace
    clear
    set obs 300000
    gen long mid = _n
    gen str10 id = "key" + string(floor(runiform() * 60000), "%07.0f")
    replace id = "" if runiform() < .01
    gen double sh = cond(runiform() < .3, ., runiform())
    gen str3 shs = cond(runiform() < .3, "", "m" + string(mod(_n, 97)))
    gen int mx = mod(_n, 30000)
    save `m', replace
}
_cjx_parity, master(`m') using(`u') keys(id) opts(unmatched(both) update) ///
    testname("large string-key join, 300k x 100k, update")
quietly {
    use `m', clear
    cjoinby id using `u', unmatched(both) update replace threads(1)
    save `ref', replace
    use `m', clear
    cjoinby id using `u', unmatched(both) update replace threads(16)
    capture cf _all using `ref', all
    local rc1 = _rc
    * the using data in another order give the same pairings
    use `u', clear
    gen double r = runiform()
    sort r
    drop r
    save `out', replace
    use `m', clear
    cjoinby id using `out', unmatched(both) update replace
    sort mid uid
    save `out', replace
    use `ref', clear
    sort mid uid
    capture cf _all using `out', all
    local rc2 = _rc
}
if !`rc1' test_pass "large join, threads(1) and threads(16) identical row for row"
else test_fail "large join, threads(1) and threads(16) identical row for row" "cf r(`rc1')"
if !`rc2' test_pass "large join, order of the using data does not change the pairings"
else test_fail "large join, order of the using data does not change the pairings" "cf r(`rc2')"
quietly {
    use `m', clear
    gen double r = runiform()
    sort r
    drop r
    save `ref', replace
}
_cjx_parity, master(`ref') using(`u') keys(id) opts(unmatched(master)) ///
    testname("large join, shuffled master, unmatched(master)")
quietly {
    clear
    set seed 17013
    set obs 1000000
    gen long id = floor(runiform() * 2000000)
    gen double mx = runiform()
    save `m', replace
    clear
    set obs 500000
    gen long id = _n * 4
    gen float ux = _n
    save `u', replace
}
_cjx_parity, master(`m') using(`u') keys(id) opts(unmatched(master)) ///
    testname("1M-observation lookup join, unmatched(master)")

print_summary "cjoinby_extra"
global CTOOLS_COMPONENT_COMPLETE "cjoinby_extra"
    global CTOOLS_CASE_COMPLETE "cjoinby_extra"
    exit
}

di as error "Unknown cjoinby case: `__case'"
exit 198
