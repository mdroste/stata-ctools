gettoken __case 0 : 0
if "`__case'" == "" local __case cdistinct
if "`__case'" == "cdistinct" {
do "validation/support/validate_setup.do"

* The reference is distinct 1.4.3 (Stata Journal 23(4), package dm0042_5).
* SSC distributes 1.2.1, whose defaults and table are different, so the suite
* stops when any other version comes first on the adopath.
capture findfile distinct.ado
if _rc {
    di as error "cdistinct suite: distinct is not installed. It needs version 1.4.3:"
    di as error "    net install dm0042_5, from(http://www.stata-journal.com/software/sj23-4)"
    exit 601
}
local __reference `"`r(fn)'"'
tempname __handle
file open `__handle' using `"`__reference'"', read text
file read `__handle' __header
file close `__handle'
if strpos(`"`__header'"', "*! 1.4.3 NJC 3 September 2023") != 1 {
    di as error `"cdistinct suite: the reference must be distinct 1.4.3, but `__reference' begins"'
    di as error `"    `__header'"'
    di as error "SSC has 1.2.1. Version 1.4.3 is in the Stata Journal archive:"
    di as error "    net install dm0042_5, from(http://www.stata-journal.com/software/sj23-4) replace"
    exit 9
}

* _cd_output runs a command and returns what it printed (r(text), lines
* joined by ;), its return code and its stored results. The timing lines of
* verbose are counted (r(timing)) and left out of the text; r(late) counts
* those printed after the table began.
capture program drop _cd_output
program define _cd_output, rclass
    _ct_loud _cd_output
    tempfile output
    tempname handle
    capture log close _cdout
    quietly log using `"`output'"', text replace name(_cdout)
    capture noisily `0'
    local rc = _rc
    local n .
    local d .
    if !`rc' {
        local n = r(N)
        local d = r(ndistinct)
    }
    quietly log close _cdout
    local text
    local timing 0
    local late 0
    local table 0
    file open `handle' using `"`output'"', read text
    file read `handle' line
    while !r(eof) {
        if substr(`"`macval(line)'"', 1, 10) == "cdistinct:" {
            local ++timing
            local late = `late' + `table'
        }
        else if substr(`"`macval(line)'"', 1, 5) != "r; t=" {
            if substr(`"`macval(line)'"', 1, 1) == "-" local table 1
            local text `"`macval(text)'`macval(line)';"'
        }
        file read `handle' line
    }
    file close `handle'
    return local text `"`macval(text)'"'
    return scalar rc = `rc'
    return scalar N = `n'
    return scalar ndistinct = `d'
    return scalar timing = `timing'
    return scalar late = `late'
end

* _cd_same passes when the reference command and the ctools command print
* the same text, end with the same return code, store the same r(N) and
* r(ndistinct), and leave the data and their sort order as they were.
capture program drop _cd_same
program define _cd_same
    args name reference candidate
    _ct_loud _cd_same
    tempfile before
    quietly save `"`before'"', emptyok
    local sorted : sortedby
    _cd_output `reference'
    local want_text `"`r(text)'"'
    local want_rc = r(rc)
    local want_n = r(N)
    local want_d = r(ndistinct)
    _cd_output `candidate'
    local got_text `"`r(text)'"'
    local got_rc = r(rc)
    local got_n = r(N)
    local got_d = r(ndistinct)
    local late = r(late)
    capture cf _all using `"`before'"'
    local same_data = !_rc
    local now : sortedby
    if "`now'" != "`sorted'" local same_data 0
    local why
    if `"`macval(want_text)'"' != `"`macval(got_text)'"' local why `why' output
    if `want_rc' != `got_rc' local why `why' rc(`want_rc' vs `got_rc')
    if `want_n' != `got_n' | `want_d' != `got_d' local why `why' r()
    if !`same_data' local why `why' data
    if `late' local why `why' timing
    if "`why'" == "" test_pass "`name'"
    else {
        di as error `"`name': differs in `why'"'
        di as error `"  want: `macval(want_text)'"'
        di as error `"  got:  `macval(got_text)'"'
        test_fail "`name'" "differs in `why'"
    }
end

* _cd_rows returns the row names of a table captured by _cd_output, in the
* order printed.
capture program drop _cd_rows
program define _cd_rows, rclass
    args text
    local names
    while `"`macval(text)'"' != "" {
        gettoken line text : text, parse(";")
        if `"`macval(line)'"' == ";" continue
        local bar = strpos(`"`macval(line)'"', "|")
        if `bar' {
            local name = strtrim(substr(`"`macval(line)'"', 1, `bar' - 1))
            if "`name'" != "" local names `names' `name'
        }
    }
    return local names `names'
end

quietly {
clear
set obs 5000
set seed 37197
gen long order = _n
gen double x = mod(_n,13)
replace x = .a if mod(_n,17)==0
replace x = .b if mod(_n,19)==0
gen str8 s = "s"+string(mod(_n,7))
replace s = "" if mod(_n,23)==0
gen byte g = mod(_n,3)
local stopped 0
local hidden 0
local shown 0
foreach keys in "x" "s" "x s" {
    foreach missing in "" "missing" {
        foreach joint in "" "joint" {
            foreach selection in "" "if g==1" "in 17/901" "if 0" {
                foreach filter in "" "minimum(0)" "minimum(8)" "maximum(7)" "minimum(14) maximum(0)" {
                    * Freeze the row restriction before the reference sorts each
                    * variable. distinct re-evaluates in after those sorts.
                    tempvar selected
                    mark `selected' `selection'
                    capture distinct `keys' if `selected', `missing' `joint' `filter'
                    local rc = _rc
                    local n = r(N)
                    local d = r(ndistinct)
                    drop `selected'
                    capture cdistinct `keys' `selection', `missing' `joint' `filter' threads(1)
                    assert _rc == `rc'
                    if !`rc' {
                        assert r(N)==`n' & r(ndistinct)==`d'
                        if "`selection'" != "if 0" & `d' == 0 local ++hidden
                        if `d' > 0 local ++shown
                    }
                    else {
                        assert `rc' == 2000 & "`joint'" != "" & "`selection'" == "if 0"
                        local ++stopped
                    }
                    assert order==_n
                }
            }
        }
    }
}
* The loop must have met each outcome: counts returned, zeros returned because
* minimum() or maximum() hid every row, and joint stopped by no observations.
assert `shown' > 0 & `hidden' > 0 & `stopped' == 30
test_pass "distinct oracle: numeric, string, joint, missing, if and in"

cdistinct x s, missing joint threads(4)
local n=r(N)
local d=r(ndistinct)
cdistinct x s, missing joint threads(1)
assert r(N)==`n' & r(ndistinct)==`d' & `d' > 0
* Results are those of the last variable while any row is shown, and zero
* once minimum() and maximum() leave no row to show.
cdistinct x s, minimum(8)
assert r(N)==4783 & r(ndistinct)==7
cdistinct s x, minimum(8)
assert r(ndistinct)==13
cdistinct x s, minimum(14)
assert r(N)==0 & r(ndistinct)==0
cdistinct x s, joint maximum(7)
assert r(N)==0 & r(ndistinct)==0
test_pass "minimum and maximum hide rows; results are zero only when no row is shown; thread parity"

sort g order
noisily _cd_same "by prefix" "by g: distinct x s, sort(distinct descending)" ///
    "by g: cdistinct x s, sort(distinct descending)"
noisily _cd_same "by prefix, joint" "by g: distinct x s, joint missing" "by g: cdistinct x s, joint missing"
noisily _cd_same "by prefix, joint, a group without observations" ///
    "by g: distinct x s if g != 1, joint" "by g: cdistinct x s if g != 1, joint"
sort order
cdistinct
assert r(N)==_N
capture cdistinct x, threads(-1)
assert _rc==198
capture cdistinct x, threads(257)
assert _rc==198
clear
set obs 0
gen x=.
* With no observations in memory distinct stops with r(111), naming a temporary
* variable that by: never created. cdistinct does what both do after if 0.
capture distinct x
assert _rc==111
noisily _cd_output cdistinct x
assert r(rc)==0 & r(N)==0 & r(ndistinct)==0
assert `"`r(text)'"' == "note: no distinct values satisfy specification;"
noisily _cd_output cdistinct x, minimum(0)
assert r(rc)==0 & r(N)==0 & r(ndistinct)==0
_cd_rows `"`r(text)'"'
assert "`r(names)'" == "x"
noisily _cd_same "empty data, joint" "distinct x, joint" "cdistinct x, joint"
capture cdistinct x, joint
assert _rc==2000
test_pass "by prefix, default varlist, validation and empty data"

* Compare everything printed with distinct: the table, the note and the
* messages. Every variable has its own number of distinct values and of
* nonmissing values, so that no sort() below meets a tie.
clear
set obs 60
gen byte short = mod(_n,3)
gen byte very_long_variable_name = mod(_n,7) if _n > 1
gen byte Zeta = mod(_n,5) if _n > 2
gen int _under = mod(_n,11) if _n > 3
gen byte alpha9 = mod(_n,2) if _n > 4
gen str4 name_s = "n"+string(mod(_n,4)) if _n > 5
local all short very_long_variable_name Zeta _under alpha9 name_s
foreach options in "" "abbrev(5)" "abbrev(32)" "abbrev(3)" "abbrev(10)" "abbrev(-2)" "abbrev(40)" ///
    "missing" "minimum(4)" "maximum(4)" "minimum(0)" "minimum(7) maximum(3)" "maximum(0)" ///
    "minimum(100)" "maximum(-1)" "maximum(-3)" "width(15)" "width(3)" "width(40)" "width(0)" ///
    "missing width(12) abbrev(8) minimum(3) maximum(6)" {
    noisily _cd_same "table: `options'" "distinct `all', `options'" "cdistinct `all', `options'"
}
test_pass "table layout, abbrev(), width(), minimum(), maximum() and the note match distinct"

foreach options in "sort(alpha)" "sort(alpha descending)" "sort(a de)" "sort(distinct)" ///
    "sort(distinct descending)" "sort(di desc)" "sort(total)" "sort(t des)" "sort(descending total)" ///
    "sort" "sort sort(descending)" "sort sort(distinct)" "sort(DIST DESC)" "sort(distinct) minimum(4)" ///
    "sort(total descending) maximum(5) abbrev(6)" "sort2(alpha)" {
    noisily _cd_same "sorted table: `options'" "distinct `all', `options'" "cdistinct `all', `options'"
}
foreach options in "sort(d)" "sort(x)" "sort(alpha total)" "sort(descending)" "sort(alphas)" ///
    "sort sort(total)" "sort(distinct alpha descending)" "sort(de d)" {
    noisily _cd_same "sort() refused: `options'" "distinct `all', `options'" "cdistinct `all', `options'"
    capture cdistinct `all', `options'
    assert _rc == 198
}
test_pass "sort() orders rows, and its errors, as distinct does"

foreach options in "joint" "joint missing" "joint minimum(100)" "joint maximum(3)" "joint sort(x)" ///
    "joint sort(descending)" "joint abbrev(5) width(20)" "joint minimum(60) maximum(1)" {
    noisily _cd_same "joint table: `options'" "distinct `all', `options'" "cdistinct `all', `options'"
}
noisily _cd_same "no observations" "distinct `all' if 0" "cdistinct `all' if 0"
noisily _cd_same "no observations, minimum(0)" "distinct `all' if 0, minimum(0) sort(alpha)" ///
    "cdistinct `all' if 0, minimum(0) sort(alpha)"
noisily _cd_same "no observations, joint" "distinct `all' if 0, joint" "cdistinct `all' if 0, joint"
capture cdistinct `all' if 0, joint
assert _rc == 2000
noisily _cd_same "if" "distinct `all' if short == 1" "cdistinct `all' if short == 1"
noisily _cd_same "one variable" "distinct Zeta" "cdistinct Zeta"
noisily _cd_same "repeated variable" "distinct Zeta short Zeta" "cdistinct Zeta short Zeta"
foreach options in "abbrev(x)" "minimum(1.5)" "width(a)" "nosuchoption" {
    noisily _cd_same "syntax error: `options'" "distinct `all', `options'" "cdistinct `all', `options'"
}
test_pass "joint table, no observations and syntax errors match distinct"

* verbose prints one timing line for each count, above the table.
noisily _cd_output distinct `all', sort(distinct)
local want `"`r(text)'"'
noisily _cd_output cdistinct `all', sort(distinct) verbose threads(2)
assert `"`r(text)'"' == `"`want'"' & r(timing) == 6 & r(late) == 0
noisily _cd_output distinct `all', joint
local want `"`r(text)'"'
noisily _cd_output cdistinct `all', joint verbose
assert `"`r(text)'"' == `"`want'"' & r(timing) == 1 & r(late) == 0
* distinct leaves its table behind as Mata globals; cdistinct creates none.
capture mata: mata drop mout
capture mata: mata drop ndistinct
capture mata: mata drop ntotal
mata: st_local("before", strofreal(rows(direxternal("*"))))
cdistinct `all', sort(total descending)
cdistinct `all', joint
mata: st_local("after", strofreal(rows(direxternal("*"))))
assert `before' == `after'
test_pass "verbose timing placement; no Mata globals are created"

* Rows that tie on the sort key. distinct orders them at random (Mata's
* order() breaks ties with the sort generator), so only the rows as a set can
* be compared with it; cdistinct keeps them in the order of varlist.
clear
set obs 24
gen byte a = mod(_n,2)
gen byte b = mod(_n,3)
gen byte c = mod(_n,2) + 5
gen byte d = mod(_n,3) + 5
gen byte e = mod(_n,4)
foreach options in "sort(distinct)" "sort(distinct descending)" "sort(total)" "sort(total descending)" {
    noisily _cd_output distinct a b c d e, `options'
    local want `"`r(text)'"'
    local want_n = r(N)
    local want_d = r(ndistinct)
    noisily _cd_output cdistinct a b c d e, `options'
    local got `"`r(text)'"'
    assert r(N) == `want_n' & r(ndistinct) == `want_d'
    mata: st_local("same", strofreal(sort(tokens(st_local("want"), ";")', 1) == sort(tokens(st_local("got"), ";")', 1)))
    assert `same'
    _cd_rows `"`got'"'
    if "`options'" == "sort(distinct)" assert "`r(names)'" == "a c b d e"
    else if "`options'" == "sort(distinct descending)" assert "`r(names)'" == "e b d a c"
    else assert "`r(names)'" == "a b c d e"
    * distinct keeps the keys in order too, whatever it does within a tie
    _cd_rows `"`want'"'
    if "`options'" == "sort(distinct)" assert inlist(word("`r(names)'", 1), "a", "c") & word("`r(names)'", 5) == "e"
    if "`options'" == "sort(distinct descending)" assert word("`r(names)'", 1) == "e"
}
noisily _cd_output cdistinct b a b e a, sort(alpha)
_cd_rows `"`r(text)'"'
assert "`r(names)'" == "a a b b e"
noisily _cd_output cdistinct b a b e a, sort(alpha descending)
_cd_rows `"`r(text)'"'
assert "`r(names)'" == "e b b a a"
test_pass "tied sort keys keep the order of varlist"

* Names beyond ASCII, and strL values the plugin cannot read directly.
clear
set obs 40
local longname = ustrunescape("gr\u00f6\u00dfenklasse_der_beobachtungen")
local cjkname = ustrunescape("\u5909\u6570\u540d\u306e\u30c6\u30b9\u30c8")
gen byte `longname' = mod(_n,8)
gen byte `cjkname' = mod(_n,5) if _n > 2
* Distinct counts 8, 5, 6 and 3, or 8, 6, 7 and 4 with missing: no ties.
gen strL big = "v" + string(mod(_n,6)) + 3000*"x" if mod(_n,9)
gen strL binary = "a" + char(0) + string(mod(_n,3)) if mod(_n,7)
foreach options in "" "abbrev(8)" "sort(alpha)" "missing sort(distinct)" "joint" "joint missing" {
    noisily _cd_same "unicode names and strL: `options'" "distinct, `options'" "cdistinct, `options'"
}
test_pass "unicode variable names and strL values"

* Above ten million observations the count columns widen with _N.
clear
set obs 10000001
gen byte x = mod(_n,3)
noisily _cd_same "more than ten million observations" "distinct x" "cdistinct x"
noisily _cd_same "more than ten million observations, width()" "distinct x, width(11) joint" "cdistinct x, width(11) joint"
test_pass "column width above ten million observations"

* Exercise hash growth, the sorting fallback, and the ordered-input path.
clear
set obs 10000
gen long order = _n
gen double x = mod(_n*7919,10000)
gen str24 s = "key"+string(x,"%012.0f")
foreach keys in "x" "s" "x s" {
    foreach cardinality in 10000 3000 100 {
        replace x = mod(order*7919,`cardinality')
        replace s = "key"+string(x,"%012.0f")
        foreach threads in 1 4 {
            cdistinct `keys', joint threads(`threads')
            assert r(N)==10000 & r(ndistinct)==`cardinality'
            assert order==_n
        }
    }
}
replace x = cond(_n<5000,0,.a)
replace s = cond(_n<5000,"", "a")
cdistinct x s, joint missing
assert r(N)==10000 & r(ndistinct)==2
capture cdistinct x s, joint
assert _rc==2000
cdistinct x s
assert r(N)==5001 & r(ndistinct)==1
test_pass "hash growth, high-cardinality fallback, and sorted missing keys"
}
noisily print_summary "cdistinct"
global CTOOLS_CASE_COMPLETE "cdistinct"
}
