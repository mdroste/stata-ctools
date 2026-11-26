* creshape: named cases share this file, each do-call has its own local scope.
gettoken __case 0 : 0
macro shift
if "`__case'" == "" local __case "creshape"

if "`__case'" == "creshape" {
* Reference comparisons with reshape (StataNow's default, favor(speed)): the
* observations in order, and everything a dataset carries besides them.
do "validation/support/validate_setup.do"
do "validation/support/rollback_helpers.do"

* Everything but the values, in global CRS_SIG: counts, sort order, data
* label, each variable's type, format, labels and characteristics, the
* dataset's characteristics and the attached value labels.
capture program drop _crs_sig
program define _crs_sig
    local sig "N=`=_N' k=`=c(k)' sort=`: sortedby' label=`: data label'"
    foreach v of varlist * {
        local chars : char `v'[]
        local chars : list sort chars
        local cs
        foreach c of local chars {
            local cs `"`cs' `c'=<`: char `v'[`c']'>"'
        }
        local sig `"`sig' | `v' `: type `v'' `: format `v'' [`: value label `v''] {`: variable label `v''}`cs'"'
    }
    local chars : char _dta[]
    local chars : list sort chars
    foreach c of local chars {
        * favor() is reshape's method; creshape records it only when given
        if "`c'" == "ReS_favor" continue
        local sig `"`sig' ; `c'=<`: char _dta[`c']'>"'
    }
    foreach v of varlist * {
        local l : value label `v'
        if "`l'" != "" {
            capture label list `l'
            if _rc == 0 local sig `"`sig' L[`l':`r(k)':`r(min)':`r(max)']"'
            else local sig `"`sig' L[`l':undefined]"'
        }
    }
    global CRS_SIG `"`sig'"'
end

* _crs_cmp "name" <arguments of reshape>: run reshape and creshape on the
* data in memory. They must return the same code; after success the data and
* their description must be equal, after a refusal creshape must leave the
* data as they were. The data in memory end as creshape left them.
capture program drop _crs_cmp
program define _crs_cmp
    gettoken name 0 : 0
    tempfile in ref
    local N0 = _N
    local sort0 : sortedby
    quietly save "`in'", emptyok orphans
    capture reshape `0'
    local rc1 = _rc
    _crs_sig
    local sig1 `"$CRS_SIG"'
    quietly save "`ref'", emptyok orphans
    quietly use "`in'", clear
    capture creshape `0'
    local rc2 = _rc
    if `rc1' != `rc2' {
        test_fail "`name'" "reshape r(`rc1'), creshape r(`rc2')"
        exit
    }
    if `rc1' {
        capture cf _all using "`in'"
        if _rc | _N != `N0' | "`: sortedby'" != "`sort0'" {
            test_fail "`name'" "data changed after r(`rc1')"
            exit
        }
        test_pass "`name'"
        exit
    }
    _crs_sig
    mata: st_local("same", strofreal(st_global("CRS_SIG") == st_local("sig1")))
    if !`same' {
        test_fail "`name'" "description differs from reshape's"
        exit
    }
    capture cf _all using "`ref'"
    if _rc {
        test_fail "`name'" "values differ from reshape's"
        exit
    }
    if c(changed) != 1 {
        test_fail "`name'" "c(changed) is `c(changed)'"
        exit
    }
    test_pass "`name'"
end

* A small wide dataset with labels, formats, notes and characteristics.
capture program drop _crs_wide
program define _crs_wide
    clear
    quietly set obs 3
    gen id = 4 - _n
    quietly replace id = 1 in 2
    quietly replace id = 2 in 3
    gen sex = id != 1
    gen int inc80 = 1000 * id
    gen long inc81 = 1100 * id
    gen double inc82 = 1200 * id + 0.5
    gen byte ue80 = id == 2
    gen byte ue81 = id == 1
    gen byte ue82 = id == 3
    label var inc80 "Income 1980"
    label var ue82 "Unemp 82"
    format inc80 %9.2f
    format inc82 %12.0gc
    label define yn 0 "no" 1 "yes"
    label values ue80 yn
    label define yn82 0 "NO" 1 "YES"
    label values ue82 yn82
    char inc80[mychar] "c80"
    char inc82[mychar] "c82"
    char _dta[mine] "dtachar"
    note inc81: note on inc81
    note inc80: note on inc80
    note inc82: note on inc80
    note: dataset note
    label data "My data"
    label var id "ID var"
    gen str5 name = "n" + string(id)
    gen extra = id * 2
    order id name
    sort sex
end

* A random wide dataset: mixed types, missing values, stubs without some
* levels, random order. The j values are left in global CRS_JLIST.
capture program drop _crs_rand
program define _crs_rand
    args seed nid J strj
    clear
    set seed `seed'
    quietly set obs `nid'
    gen long id = _n * 3 - 7
    gen str6 grp = word("a bb ccc dd e", 1 + mod(_n, 5))
    if runiform() < 0.3 quietly replace id = . in 1
    gen double c1 = runiform()
    gen byte c2 = mod(_n, 4)
    gen str10 c3 = "s" + string(mod(_n * 7, 13))
    if runiform() < 0.5 quietly replace c3 = "" if mod(_n, 5) == 0
    local jlist
    forvalues j = 1/`J' {
        if `strj' local jj : word `j' of q w e r t y u i o p aa bb cc dd ee ff gg hh
        else local jj = `j' * `j' + cond(runiform() < 0.3, 100, 0)
        local jlist `jlist' `jj'
        if runiform() < 0.85 {
            local t : word `=1 + int(5 * runiform())' of byte int long float double
            gen `t' x`jj' = floor(100 * runiform())
            if runiform() < 0.5 quietly replace x`jj' = . if runiform() < 0.2
            if runiform() < 0.3 quietly replace x`jj' = .b if runiform() < 0.1
            if runiform() < 0.3 format x`jj' %9.2f
            if runiform() < 0.3 label var x`jj' "x at `jj'"
            if runiform() < 0.3 note x`jj': note `jj'
            if runiform() < 0.3 char x`jj'[ch] "v`jj'"
        }
        if runiform() < 0.85 gen double y`jj' = rnormal()
        if runiform() < 0.85 {
            local w = 1 + int(12 * runiform())
            gen str`w' z`jj' = substr("abcdefghijklmnop", 1, int(`w' * runiform()) + 1)
            if runiform() < 0.5 quietly replace z`jj' = "" if runiform() < 0.2
        }
    }
    gen u = runiform()
    sort u
    drop u
    if runiform() < 0.5 {
        label define lab1 0 "zero" 1 "one" 2 "two"
        label values c2 lab1
    }
    if runiform() < 0.5 label data "random data"
    if runiform() < 0.5 note: a data note
    global CRS_JLIST `jlist'
end

* Long panels with a numeric, string or double id (with missing ids), some
* constant variables, an unbalanced design and a random order.
capture program drop _crs_long
program define _crs_long
    args seed n J kind
    clear
    set seed `seed'
    quietly set obs `n'
    gen long id = ceil(_n / `J')
    gen int t = mod(_n - 1, `J') + 1
    if "`kind'" == "strid" {
        gen str8 sid = "k" + string(mod(id * 7919, 100003))
        drop id
        rename sid id
    }
    if "`kind'" == "dblid" {
        gen double did = id / 7 + cond(mod(id, 11) == 0, ., 0)
        quietly replace did = .a if mod(id, 13) == 0 & missing(did)
        drop id
        rename did id
        quietly bysort id t: keep if _n == 1
    }
    gen double a = rnormal()
    gen float b = runiform()
    gen long c = floor(1e6 * runiform())
    gen str7 s = substr("zyxwvuts", 1, 1 + mod(_n, 7))
    bysort id (t): gen double k2 = a[1]
    by id: gen str5 k3 = s[1]
    quietly drop if runiform() < 0.15
    gen double u = runiform()
    sort u
    drop u
end

* ---------------------------------------------------------------------------
* Round trips on labelled data
* ---------------------------------------------------------------------------
_crs_wide
_crs_cmp "long and wide round trips equal reshape" long inc ue, i(id) j(year)
_crs_cmp "wide, settings from the data" wide
_crs_cmp "long again: notes and characteristics merged once" long
_crs_wide
_crs_cmp "listed j values, one without variables" long inc ue, i(id) j(year 80 82 85)
_crs_wide
_crs_cmp "string j from the variable names" long inc ue, i(id) j(year) string
_crs_cmp "string j back to wide" wide
_crs_wide
_crs_cmp "two i variables" long inc ue, i(sex id) j(year)
_crs_wide
_crs_cmp "one stub of two" long inc, i(id) j(year)
_crs_wide
_crs_cmp "default name of j" long inc ue, i(id)

* a value label and variable labels of the long form survive wide and back
clear
quietly set obs 30
gen id = ceil(_n / 3)
gen byte wave = mod(_n, 3) + 1
label define wv 1 "first" 2 `"the "second""' 3 "third"
label values wave wv
label variable wave "survey wave"
gen double income = _n * 1.5
label variable income "household income"
gen averyveryverylongvariablename = _n
label variable averyveryverylongvariablename "a long name"
_crs_cmp "labels set aside by wide" wide income averyveryverylongvariablename, i(id) j(wave)
_crs_cmp "labels put back by long, format widened to the labels" long

* ---------------------------------------------------------------------------
* Random data, small (serial paths)
* ---------------------------------------------------------------------------
forvalues r = 1/24 {
    local seed = 7000 + `r'
    set seed `seed'
    local nid = 1 + int(60 * runiform())
    local J = 1 + int(6 * runiform())
    local strj = runiform() < 0.3
    local sopt = cond(`strj', "string", "")
    _crs_rand `seed' `nid' `J' `strj'
    local which = mod(`r', 4)
    if `which' == 0 _crs_cmp "random `r': long i(id)" long x y z, i(id) j(t) `sopt'
    else if `which' == 1 _crs_cmp "random `r': long i(grp id)" long x y z, i(grp id) j(t) `sopt'
    else if `which' == 2 _crs_cmp "random `r': long x" long x, i(id) j(t) `sopt'
    else {
        local first : word 1 of $CRS_JLIST
        _crs_cmp "random `r': long, j listed" long x z, i(id grp) j(t `first' `=cond(`strj', "zz", "999")') `sopt'
    }
    capture confirm variable t
    if !_rc & _N {
        _crs_cmp "random `r': wide" wide
        _crs_cmp "random `r': long again" long
        set seed `seed'
        quietly drop if runiform() < 0.2
        gen u = runiform()
        sort u
        drop u
        * (with j listed, reshape reports "no data" for every value after
        * one without observations and drops their data; see below)
        if _N & `which' != 3 _crs_cmp "random `r': wide, unbalanced" wide
    }
}

* ---------------------------------------------------------------------------
* Larger data (parallel paths), sorted and unsorted, three kinds of id
* ---------------------------------------------------------------------------
foreach kind in plain strid dblid {
    _crs_long 11 120000 6 `kind'
    _crs_cmp "wide, `kind' id, random order" wide a b c s, i(id) j(t)
    _crs_cmp "long, `kind' id" long
    _crs_long 12 120000 4 `kind'
    sort id t
    _crs_cmp "wide, `kind' id, sorted" wide a b c s, i(id) j(t)
}
_crs_long 13 120000 5 plain
preserve
quietly reshape wide a b c s, i(id) j(t)
tempfile wideref
quietly save "`wideref'"
restore
quietly creshape wide a b c s, i(id) j(t) threads(1)
capture cf _all using "`wideref'"
if _rc == 0 test_pass "threads(1) equals the default"
else test_fail "threads(1) equals the default" "values differ"
quietly creshape long, threads(2)
quietly creshape wide, threads(3)
capture cf _all using "`wideref'"
if _rc == 0 test_pass "threads() with the settings from the data"
else test_fail "threads() with the settings from the data" "values differ"

* string j in the observations; listed values in and out of order
clear
set seed 5
quietly set obs 3000
gen id = ceil(_n / 3)
gen str4 w = word("pre mid post", mod(_n, 3) + 1)
gen double v = rnormal()
gen str3 q = string(mod(_n, 9))
gen u = runiform()
sort u
drop u
_crs_cmp "string j in the observations" wide v q, i(id) j(w) string
_crs_cmp "string j back to long" long
_crs_long 21 600 6 plain
_crs_cmp "wide, j listed, one value without observations" wide a b c s, i(id) j(t 2 5 9)
_crs_long 21 600 6 plain
_crs_cmp "wide, j listed out of order" wide a s, i(id) j(t 5 1 3)

* ---------------------------------------------------------------------------
* Refused requests: reshape's return code, and the data untouched
* ---------------------------------------------------------------------------
_crs_wide
_crs_cmp "no variables for the stub: r(111)" long zz, i(id) j(year)
_crs_wide
quietly replace id = 1
_crs_cmp "i does not identify the wide observations: r(9)" long inc ue, i(id) j(year)
_crs_long 22 600 4 plain
quietly replace t = 1 in 5
_crs_cmp "j repeats within i: r(9)" wide a b c s, i(id) j(t)
_crs_long 22 600 4 plain
quietly replace k2 = 5 in 7
_crs_cmp "another variable varies within i: r(9)" wide a b c s, i(id) j(t)
_crs_long 22 600 4 plain
quietly replace t = . in 3
_crs_cmp "j missing: r(498)" wide a b c s, i(id) j(t)
_crs_long 22 600 4 plain
gen a2 = 1
_crs_cmp "name of a new variable taken: r(110)" wide a b c s, i(id) j(t)
_crs_long 22 600 4 plain
_crs_cmp "stub is not a variable: r(111)" wide a zz, i(id) j(t)
_crs_long 22 600 4 plain
_crs_cmp "numeric j with string: r(109)" wide a b c s, i(id) j(t) string
_crs_long 22 600 4 plain
_crs_cmp "wide without j(): r(198)" wide a b c s, i(id)
_crs_long 22 600 4 plain
_crs_cmp "without i(): r(198)" wide a b c s, j(t)
_crs_long 22 600 4 plain
_crs_cmp "long with an existing j: r(110)" long a, i(id) j(t)
_crs_long 22 600 4 plain
quietly replace t = t + 0.5
_crs_cmp "fractional j makes invalid names: r(198)" wide a b c s, i(id) j(t)
clear
quietly set obs `=c(maxvar) + 10'
gen id = 1
gen long t = _n
gen byte a = 1
_crs_cmp "more values of j than maxvar has room for: r(134)" wide a, i(id) j(t)
clear
quietly set obs 5
gen id = _n
gen a1 = 1
gen a2 = "x"
_crs_cmp "string and numeric variables of one stub: r(198)" long a, i(id) j(t)
clear
quietly set obs 5
gen id = _n
gen a1 = 1
gen a = 2
_crs_cmp "long name taken: r(110)" long a, i(id) j(t)
clear
quietly set obs 4
gen id = ceil(_n / 2)
gen str5 t = cond(mod(_n, 2), "a", "b c")
gen a = _n
_crs_cmp "string j with a space: r(111)" wide a, i(id) j(t) string
quietly replace t = " b" if t == "b c"
_crs_cmp "string j with a leading blank: r(498)" wide a, i(id) j(t) string
quietly replace t = "" if t == " b"
_crs_cmp "string j empty: r(498)" wide a, i(id) j(t) string
_crs_cmp "string j without string: r(109)" wide a, i(id) j(t)

* messages go through the error channel: shown under quietly, as plain text
clear
quietly set obs 6
gen id = ceil(_n / 3)
gen t = mod(_n, 3) + 1
gen a = _n
gen c = _n
tempname lg
tempfile lf
quietly log using "`lf'", text name(`lg') replace
capture noisily quietly creshape wide a, i(id) j(t)
quietly log close `lg'
tempname fh
file open `fh' using "`lf'", read text
file read `fh' line
local seen_first 0
local seen_markup 0
while r(eof) == 0 {
    if strpos(`"`macval(line)'"', "variable c not constant within id") local seen_first 1
    if strpos(`"`macval(line)'"', "{bf") | strpos(`"`macval(line)'"', "{p ") local seen_markup 1
    file read `fh' line
}
file close `fh'
if `seen_first' & !`seen_markup' test_pass "refusal explained under quietly, without markup"
else test_fail "refusal explained under quietly, without markup" "first=`seen_first' markup=`seen_markup'"

* ---------------------------------------------------------------------------
* Advanced syntax: declared settings, an xi list, atwl()
* ---------------------------------------------------------------------------
clear
quietly set obs 40
gen id = _n
gen g = mod(_n, 3)
forvalues j = 1/3 {
    gen inc`j'r = _n * `j'
    gen ue`j' = _n + `j'
}
gen keepme = 3
gen dropme = 4
tempfile adv
quietly save "`adv'"
foreach cmd in reshape creshape {
    quietly use "`adv'", clear
    `cmd' clear
    `cmd' i id
    `cmd' j yr
    `cmd' xij inc@r ue, atwl(X)
    `cmd' xi keepme g
    quietly `cmd' long
    _crs_sig
    local sigl_`cmd' `"$CRS_SIG"'
    tempfile outl_`cmd'
    quietly save "`outl_`cmd''"
    quietly `cmd' wide
    _crs_sig
    local sigw_`cmd' `"$CRS_SIG"'
    tempfile outw_`cmd'
    quietly save "`outw_`cmd''"
}
quietly use "`outl_creshape'", clear
capture cf _all using "`outl_reshape'"
if _rc == 0 & `"`sigl_reshape'"' == `"`sigl_creshape'"' test_pass "declared settings, xi list in its order, atwl: long"
else test_fail "declared settings, xi list in its order, atwl: long" "differs from reshape"
quietly use "`outw_creshape'", clear
capture cf _all using "`outw_reshape'"
if _rc == 0 & `"`sigw_reshape'"' == `"`sigw_creshape'"' test_pass "declared settings, xi list, atwl: wide"
else test_fail "declared settings, xi list, atwl: wide" "differs from reshape"

* atwl() in the basic syntax (reshape takes it only in reshape xij)
quietly use "`adv'", clear
quietly creshape long inc@r ue, i(id) j(yr) atwl(X)
capture confirm variable incXr ueX
if _rc == 0 & "`: char _dta[ReS_atwl]'" == "X" test_pass "atwl() in the basic syntax"
else test_fail "atwl() in the basic syntax" "variables incXr ueX not made"
quietly reshape wide
capture confirm variable inc1r inc2r inc3r ue1 ue2 ue3
if _rc == 0 test_pass "reshape continues from creshape's settings"
else test_fail "reshape continues from creshape's settings" "r(`=_rc')"

* favor() is reshape's choice of method: recorded for reshape, not used
quietly use "`adv'", clear
quietly creshape long inc@r ue, i(id) j(yr) favor(memory)
local fav : char _dta[ReS_favor]
capture creshape wide inc@r ue, i(id) j(yr) favor(fast)
if "`fav'" == "memory" & _rc == 198 test_pass "favor() is recorded and checked"
else test_fail "favor() is recorded and checked" "char `fav', r(`=_rc')"

* query, clear and error are reshape's
quietly use "`adv'", clear
capture creshape query
local rc_q = _rc
quietly creshape long inc@r ue, i(id) j(yr)
capture creshape error
local rc_e = _rc
creshape clear
if `rc_q' == 0 & `rc_e' == 0 & "`: char _dta[ReS_i]'" == "" test_pass "query, error and clear"
else test_fail "query, error and clear" "r(`rc_q') r(`rc_e')"
capture creshape long
if _rc == 111 test_pass "no settings: reshape's r(111)"
else test_fail "no settings: reshape's r(111)" "r(`=_rc')"

* ---------------------------------------------------------------------------
* Results, timing table, nothing left behind
* ---------------------------------------------------------------------------
_crs_wide
creshape long inc ue, i(id) j(year)
if r(N) == 9 & r(k) == 7 & r(N_before) == 3 & r(k_before) == 10 & "`r(cmd)'" == "creshape" & ///
    r(time_store) < . & r(threads_max) >= 1 {
    test_pass "stored results"
}
else test_fail "stored results" "r(N)=`r(N)' r(k)=`r(k)'"
quietly log using "`lf'", text name(`lg') replace
creshape wide, verbose
quietly log close `lg'
file open `fh' using "`lf'", read text
file read `fh' line
local seen 0
while r(eof) == 0 {
    if strpos(`"`macval(line)'"', "creshape timing breakdown") local ++seen
    if strpos(`"`macval(line)'"', "Write data to Stata") local ++seen
    if strpos(`"`macval(line)'"', "Groups and cells") local ++seen
    file read `fh' line
}
file close `fh'
capture scalar list _creshape_time_load
if `seen' == 3 & _rc test_pass "verbose timing table; no scalars are left"
else test_fail "verbose timing table; no scalars are left" "lines `seen', scalar r(`=_rc')"

* ---------------------------------------------------------------------------
* A failure after the data started to change puts them back
* ---------------------------------------------------------------------------
capture program drop _crs_fixture
program define _crs_fixture
    clear
    set seed 7
    quietly set obs 500
    gen long id = _n
    gen str9 name = "n" + string(_n)
    forvalues j = 1/4 {
        gen double a`j' = rnormal()
        gen str5 s`j' = string(`j' * _n)
    }
    label define yn 0 "no" 1 "yes"
    gen byte flag = mod(_n, 2)
    label values flag yn
    note a1: a note
    char a2[c] "x"
    char _dta[d] "dd"
    label data "rollback data"
    sort name
end
foreach fault in template schema store {
    _crs_fixture
    * the settings are declared first: like reshape, creshape keeps them
    * with the data when a request fails
    creshape i id
    creshape j t
    creshape xij a s
    _ctools_rb_case, cmd(creshape long) global(CTOOLS_creshape_fault) ///
        fault(`fault') testname("long put back after a failure at `fault'") modes(default)
    quietly reshape long
    quietly drop if mod(_n, 7) == 0
    sort name t
    _ctools_rb_case, cmd(creshape wide) global(CTOOLS_creshape_fault) ///
        fault(`fault') testname("wide put back after a failure at `fault'") modes(default)
}
* a saved dataset stays unchanged and keeps its file name; the next call works
_crs_fixture
tempfile disk
quietly save "`disk'"
global CTOOLS_creshape_fault store
capture creshape long a s, i(id) j(t)
global CTOOLS_creshape_fault
if _rc == 1 & c(changed) == 0 & c(filename) == "`disk'" test_pass "c(changed) and c(filename) after a failure"
else test_fail "c(changed) and c(filename) after a failure" "r(`=_rc') changed=`c(changed)'"
quietly creshape long a s, i(id) j(t)
if _N == 2000 & c(changed) == 1 & c(filename) == "" test_pass "the next call starts clean"
else test_fail "the next call starts clean" "N=`=_N'"

* ---------------------------------------------------------------------------
* Edges
* ---------------------------------------------------------------------------
* strL variables: reshape does the work
_crs_fixture
gen strL big = name + "!"
_crs_cmp "strL data: long through reshape" long a s, i(id) j(t)
_crs_cmp "strL data: wide through reshape" wide

* no observations (reshape itself stops with a Mata error here)
_crs_fixture
quietly drop in 1/l
capture creshape long a s, i(id) j(t)
local ok = _rc == 0 & _N == 0 & c(k) == 6 & "`: sortedby'" == "id t"
capture confirm variable id t a s name flag
if `ok' & _rc == 0 test_pass "no observations: long"
else test_fail "no observations: long" "N=`=_N' k=`=c(k)'"
capture creshape wide
if _rc == 0 & _N == 0 & "`: sortedby'" == "id" test_pass "no observations: wide"
else test_fail "no observations: wide" "r(`=_rc')"

clear
quietly set obs 1
gen id = 5
gen x3 = 2
_crs_cmp "one observation, one level: long" long x, i(id) j(t)
_crs_cmp "one observation, one level: wide" wide

* many variables
clear
set seed 3
quietly set obs 200
gen long id = _n
forvalues j = 1/400 {
    gen float v`j' = runiform()
    gen byte w`j' = mod(_n + `j', 5)
}
gen u = runiform()
sort u
drop u
_crs_cmp "800 wide variables: long" long v w, i(id) j(t)
_crs_cmp "400 levels: wide" wide

* names and strings outside ASCII; i with missing values, blanks, -0, 1e300
clear
set seed 31
quietly set obs 12
gen str12 gruppe = word("Zebra apfel Äpfel zebra 10 9 a_b éclair", 1 + mod(_n, 8))
quietly replace gruppe = "" in 3
quietly replace gruppe = "a b" in 4
gen double key = cond(mod(_n, 3) == 0, ., _n / 3)
quietly replace key = .a in 6
quietly replace key = .z in 12
quietly replace key = -0 in 1
quietly replace key = 1e300 in 2
gen größe1 = _n
gen größe2 = -_n
gen str8 naïve1 = "é" + string(_n)
gen str8 naïve2 = "ü"
gen u = runiform()
sort u
drop u
_crs_cmp "names outside ASCII, odd i values: long" long größe naïve, i(gruppe key) j(jähr)
_crs_cmp "names outside ASCII, odd i values: wide" wide

* the sort marker is true: by works, and a real sort gives the same order
_crs_fixture
quietly creshape long a s, i(name id) j(t)
capture noisily {
    by name id (t): assert _n == t
    gen long before = _n
    gen double shuffle = runiform()
    sort shuffle
    sort name id t
    assert before == _n
    drop before shuffle
}
if _rc == 0 test_pass "long result is in the order its sort marker says"
else test_fail "long result is in the order its sort marker says" "r(`=_rc')"
quietly creshape wide
capture noisily {
    by name: assert _N == 1
    gen long row = _n
    sort id
    sort name id
    assert row == _n
}
if _rc == 0 test_pass "wide result is in the order its sort marker says"
else test_fail "wide result is in the order its sort marker says" "r(`=_rc')"

* ---------------------------------------------------------------------------
* Where creshape departs from reshape on purpose
* ---------------------------------------------------------------------------
* a j value beyond the range of a long (reshape stores it as missing)
clear
quietly set obs 3
gen id = _n
gen a100 = 1
gen a3000000000 = 3
quietly creshape long a, i(id) j(t)
capture assert inlist(t, 100, 3000000000)
if _rc == 0 & "`: type t'" == "double" test_pass "large j values are kept"
else test_fail "large j values are kept" "type `: type t'"

* j listed, the first value without observations: the others keep theirs
clear
input id str2 t x
1 "zz" 2
2 "zz" 3
end
quietly creshape wide x, i(id) j(t q zz) string
capture assert xzz == id + 1
if _rc == 0 test_pass "listed j: only values without observations are left out"
else test_fail "listed j: only values without observations are left out" "r(`=_rc')"

* value labels not attached to a variable are kept
_crs_wide
label define orphan 7 "seven"
quietly creshape long inc ue, i(id) j(year)
capture label list orphan
if _rc == 0 test_pass "unattached value labels are kept"
else test_fail "unattached value labels are kept" "label dropped"

print_summary "creshape"
global CTOOLS_COMPONENT_COMPLETE "creshape"
    global CTOOLS_CASE_COMPLETE "creshape"
    exit
}

di as error "Unknown creshape case: `__case'"
exit 198
