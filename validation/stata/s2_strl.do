* s2_strl: named cases share this file, each do-call has its own local scope.
gettoken __case 0 : 0
macro shift
if "`__case'" == "" local __case "s2_strl"

if "`__case'" == "s2_strl" {
do validation/support/validate_setup.do
capture confirm file "validation/support/rollback_helpers.do"
if !_rc do "validation/support/rollback_helpers.do"
else do "rollback_helpers.do"
set seed 20261104

* strL values the plugin cannot read: text over 2,045 bytes (3,000 and 10,000
* bytes, some differing only in the last byte), binary values from fileread()
* (with and without NUL bytes), plus "", short text and missing numbers.
capture program drop _s2s_data
program define _s2s_data
    syntax [, n(integer 210)]
    clear
    quietly set obs `n'
    tempfile b1 b2
    tempname fh
    file open `fh' using `b1', write binary replace
    file write `fh' %1bu (0) %1bu (65) %1bu (0) %1bu (66)
    file close `fh'
    file open `fh' using `b2', write binary replace
    file write `fh' %1bu (200) %1bu (49) %1bu (255) %1bu (32)
    file close `fh'
    generate long id = _n
    * by() key: nine values in a scrambled order of observations
    generate int k = mod(_n * 7, 9)
    generate strL g = ""
    quietly replace g = "A" * 3000 if k == 1
    quietly replace g = "A" * 2999 + "B" if k == 2
    quietly replace g = "A" * 10000 if k == 3
    quietly replace g = "A" * 9999 + "C" if k == 4
    quietly replace g = fileread("`b1'") if k == 5
    quietly replace g = fileread("`b2'") if k == 6
    quietly replace g = "x" if k == 7
    quietly replace g = "A" * 2999 if k == 8
    * a second strL key the plugin can read, and a fixed string key
    generate strL h = "h" + string(mod(_n, 4))
    generate str2 s = "s" + string(mod(_n, 3))
    * strL payload: long, binary and empty
    generate strL note = "n" * 3000 + string(_n)
    quietly replace note = fileread("`b1'") if mod(_n, 17) == 0
    quietly replace note = fileread("`b2'") if mod(_n, 19) == 0
    quietly replace note = "" if mod(_n, 23) == 0
    label variable note "payload"
    note note: a note on the payload
    generate double x = mod(_n * 13, 31)
    quietly replace x = . if mod(_n, 29) == 0
    generate double y = _n / 4 + mod(_n, 5)
    quietly replace y = . if mod(_n, 11) == 0
    drop k
end

* Stata's sort rank of a strL key as a fixed string, the key the plugin
* would group and order exactly as it does v; with empty, "" stays ""
capture program drop _s2s_rankkey
program define _s2s_rankkey
    args v newvar empty
    tempvar r
    quietly egen long `r' = group(`v'), missing
    quietly generate str6 `newvar' = string(`r', "%06.0f")
    if "`empty'" != "" quietly replace `newvar' = "" if `v' == ""
end

* The error lines of two text logs agree (commands, the log header and
* rangejoin's note on rangestat's version left out)
capture mata: mata drop _s2s_loglines()
capture mata: mata drop _s2s_errlines()
mata:
string colvector _s2s_errlines(string scalar f)
{
    string colvector l
    real scalar i
    real colvector keep

    l = cat(f)
    keep = J(rows(l), 1, 1)
    for (i = 1; i <= rows(l); i++) {
        if (strtrim(l[i]) == "" | substr(l[i], 1, 2) == ". " | substr(l[i], 1, 2) == "> " |
            strpos(l[i], "(using rangestat") | strpos(l[i], "name:") | strpos(l[i], "log:") |
            strpos(l[i], "log type:") | strpos(l[i], "opened on:") | strpos(l[i], "closed on:") |
            regexm(l[i], "^-+$")) keep[i] = 0
    }
    return(select(l, keep))
}

void _s2s_loglines(string scalar f1, string scalar f2)
{
    string colvector a, b
    a = _s2s_errlines(f1)
    b = _s2s_errlines(f2)
    if (rows(a) == 0 | rows(a) != rows(b)) exit(error(9))
    if (any(a :!= b)) exit(error(9))
}
end

_s2s_data
tempfile base
quietly save `base'
quietly count if strlen(g) > 2045
assert r(N) > 0

* The plugin's check: 1 for a strL variable with a value it cannot load
capture noisily {
    use `base', clear
    _ctools_load
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    local ctools_strl_bad
    plugin call ctools_plugin id g h s note x, "strlcheck"
    assert "`ctools_strl_bad'" == "0 1 0 0 1 0"
    * only the in-range counts: observations 7 to 8 hold short text in g
    quietly replace g = "y" in 7/8
    plugin call ctools_plugin g in 7/8, "strlcheck"
    assert "`ctools_strl_bad'" == "0"
    _ctools_strlkey id g h s note, check
    assert "`__ctools_strl_bad'" == "g note"
}
if _rc test_fail "strlcheck finds the strL values the plugin cannot load" "r(`=_rc')"
else test_pass "strlcheck finds the strL values the plugin cannot load"

* ---------------------------------------------------------------------------
* SORT-03: csort with strL variables (keys and payload) is sort, stable
* ---------------------------------------------------------------------------
foreach keys in "g" "g x" "x" "s x" "note" "h g s" "y g" {
    capture noisily {
        use `base', clear
        sort `keys', stable
        tempfile expected
        quietly save `expected'
        local marker : sortedby
        use `base', clear
        csort `keys'
        assert r(N) == _N & "`r(cmd)'" == "csort" & "`r(sortvars)'" == "`keys'"
        cf _all using `expected'
        assert "`: sortedby'" == "`marker'"
        use `base', clear
        csort `keys', nosortedby
        cf _all using `expected'
        assert "`: sortedby'" == ""
    }
    if _rc test_fail "SORT-03 csort `keys' with strL variables" "r(`=_rc')"
    else test_pass "SORT-03 csort `keys' with strL variables"
}
* Stable: ties keep their order (id ascending within key), e(sample) moves
capture noisily {
    use `base', clear
    quietly regress y x if mod(id, 3)
    csort g
    assert id > id[_n - 1] if g == g[_n - 1] & _n > 1
    generate byte es = e(sample)
    tempfile got
    quietly save `got'
    use `base', clear
    quietly regress y x if mod(id, 3)
    sort g, stable
    generate byte es = e(sample)
    cf _all using `got'
    * already sorted: nothing to do, and still sorted
    csort g
    cf _all using `got'
}
if _rc test_fail "SORT-03 csort strL: stable ties and e(sample)" "r(`=_rc')"
else test_pass "SORT-03 csort strL: stable ties and e(sample)"

* ---------------------------------------------------------------------------
* DIST-01: cdistinct counts long and binary strL values as distinct does
* ---------------------------------------------------------------------------
foreach spec in "g" "g, missing" "g h s" "g x, joint" "g s, joint missing" ///
    "g note, missing" "g if x > 10" "g in 20/150, missing" "note" {
    capture noisily {
        use `base', clear
        distinct `spec'
        local rn = r(N)
        local rd = r(ndistinct)
        cdistinct `spec'
        assert r(N) == `rn' & r(ndistinct) == `rd'
    }
    if _rc test_fail "DIST-01 cdistinct `spec'" "r(`=_rc')"
    else test_pass "DIST-01 cdistinct `spec'"
}
* values that differ only in their last byte count separately
capture noisily {
    clear
    quietly set obs 4
    generate strL t = "A" * 5000
    quietly replace t = "A" * 4999 + "B" in 2
    quietly replace t = "A" * 4999 + "B" in 3
    cdistinct t
    assert r(N) == 4 & r(ndistinct) == 2
    distinct t
    assert r(ndistinct) == 2
}
if _rc test_fail "DIST-01 cdistinct long values differing at the end" "r(`=_rc')"
else test_pass "DIST-01 cdistinct long values differing at the end"
* by: (recall) with a long strL variable: r() of the last group
capture noisily {
    use `base', clear
    distinct g if s == "s2", missing
    local rn = r(N)
    local rd = r(ndistinct)
    sort s
    by s: cdistinct g, missing
    assert r(N) == `rn' & r(ndistinct) == `rd'
}
if _rc test_fail "DIST-01 by: cdistinct with long strL" "r(`=_rc')"
else test_pass "DIST-01 by: cdistinct with long strL"

* ---------------------------------------------------------------------------
* IPOL-03: cinterpolate by() with long and binary strL keys is ipolate
* ---------------------------------------------------------------------------
foreach opts in "" "epolate" {
    foreach by in "g" "g s" "h g" {
        capture noisily {
            use `base', clear
            ipolate y x, by(`by') generate(expected) `opts'
            cinterpolate y x, by(`by') generate(got) `opts'
            assert expected == got
            assert "`: type got'" == "`: type expected'"
            * by: prefix
            use `base', clear
            sort `by'
            by `by': ipolate y x, generate(expected) `opts'
            by `by': cinterpolate y x, generate(got) `opts'
            assert expected == got
            * if/in
            use `base', clear
            ipolate y x if id > 30 in 1/180, by(`by') generate(expected) `opts'
            cinterpolate y x if id > 30 in 1/180, by(`by') generate(got) `opts'
            assert expected == got
        }
        if _rc test_fail "IPOL-03 cinterpolate by(`by') `opts'" "r(`=_rc')"
        else test_pass "IPOL-03 cinterpolate by(`by') `opts'"
    }
}

* ---------------------------------------------------------------------------
* RSTAT-07: crangestat by() with long and binary strL keys is rangestat
* ---------------------------------------------------------------------------
foreach by in "g" "g s" "h g" {
    foreach xs in "" "excludeself" {
        capture noisily {
            use `base', clear
            rangestat (count) y (sum) y (mean) y (min) y (max) y (sd) y, ///
                interval(x -3 3) by(`by') `xs'
            tempfile expected
            quietly save `expected'
            use `base', clear
            crangestat (count) y (sum) y (mean) y (min) y (max) y (sd) y, ///
                interval(x -3 3) by(`by') `xs'
            foreach v in y_count y_sum y_mean y_min y_max y_sd {
                rename `v' got_`v'
            }
            quietly merge 1:1 id using `expected', keepusing(y_*) nogenerate
            foreach v in y_count y_sum y_mean y_min y_max {
                assert got_`v' == `v'
            }
            assert reldif(got_y_sd, y_sd) < 1e-9 | (got_y_sd == . & y_sd == .)
        }
        if _rc test_fail "RSTAT-07 crangestat by(`by') `xs'" "r(`=_rc')"
        else test_pass "RSTAT-07 crangestat by(`by') `xs'"
    }
}
* the flexible statistics (Mata) group by the strL variable itself
capture noisily {
    use `base', clear
    rangestat (cov) x y, interval(x -5 5) by(g)
    tempfile expected
    quietly save `expected'
    use `base', clear
    crangestat (cov) x y (sum) y, interval(x -5 5) by(g)
    drop y_sum
    cf _all using `expected'
}
if _rc test_fail "RSTAT-07 crangestat (cov) with (sum) by(long strL)" "r(`=_rc')"
else test_pass "RSTAT-07 crangestat (cov) with (sum) by(long strL)"
* rollback: a failure after the outputs exist drops them (stand-ins too)
use `base', clear
_ctools_rb_case, cmd(crangestat (sum) y, interval(x -2 2) by(g)) ///
    global(CTOOLS_crangestat_fault) fault(outputs) testname("RSTAT-07 rollback with strL by()")
capture noisily {
    use `base', clear
    local k0 = c(k)
    crangestat (sum) y, interval(x -2 2) by(g)
    assert c(k) == `k0' + 1
    confirm variable y_sum
}
if _rc test_fail "RSTAT-07 crangestat leaves only its result variable" "r(`=_rc')"
else test_pass "RSTAT-07 crangestat leaves only its result variable"

* ---------------------------------------------------------------------------
* SAMP-06: csample by() with long and binary strL keys
* ---------------------------------------------------------------------------
* The draw matches the draw for a fixed-string key with the same groups in
* the same order (the plugin's own path), and the kept count per group
* matches sample's.
foreach spec in "50|" "33.3|" "4|count" "1|count" {
    gettoken pct rest : spec, parse("|")
    local cnt = cond(strpos("`spec'", "count"), "count(`pct')", "")
    local arg = cond("`cnt'" != "", "", "`pct'")
    local native = cond("`cnt'" != "", "`pct', count", "`pct',")
    capture noisily {
        use `base', clear
        _s2s_rankkey g rk
        set seed 99
        csample `arg', `cnt' by(rk)
        keep id
        tempfile expected
        quietly save `expected'
        use `base', clear
        set seed 99
        csample `arg', `cnt' by(g)
        keep id
        cf _all using `expected'
        * per-group counts: sample
        use `base', clear
        _s2s_rankkey g rk
        set seed 5
        sample `native' by(g)
        contract rk, freq(n_native)
        tempfile native
        quietly save `native'
        use `base', clear
        _s2s_rankkey g rk
        set seed 5
        csample `arg', `cnt' by(g)
        contract rk, freq(n_native)
        cf _all using `native'
    }
    if _rc test_fail "SAMP-06 csample `spec' by(long strL)" "r(`=_rc')"
    else test_pass "SAMP-06 csample `spec' by(long strL)"
}
capture noisily {
    use `base', clear
    _s2s_rankkey g rk
    set seed 7
    csample 40 if x > 5, by(rk s)
    keep id
    tempfile expected
    quietly save `expected'
    use `base', clear
    set seed 7
    csample 40 if x > 5, by(g s)
    keep id
    cf _all using `expected'
    * by: prefix
    use `base', clear
    _s2s_rankkey g rk
    sort rk, stable
    set seed 8
    by rk: csample 50
    keep id
    quietly save `expected', replace
    use `base', clear
    sort g, stable
    set seed 8
    by g: csample 50
    keep id
    cf _all using `expected'
}
if _rc test_fail "SAMP-06 csample if, by(strL fixed), by: prefix" "r(`=_rc')"
else test_pass "SAMP-06 csample if, by(strL fixed), by: prefix"

* ---------------------------------------------------------------------------
* BSAMP-01: cbsample strata()/cluster() with long and binary strL keys
* ---------------------------------------------------------------------------
foreach opts in "strata(g)" "cluster(g)" "strata(s) cluster(g)" "strata(g) cluster(note)" ///
    "cluster(g) idcluster(newid)" "strata(g) weight(w)" "3, strata(s) cluster(g)" {
    capture noisily {
        local n
        local o `opts'
        if strpos("`opts'", ",") {
            gettoken n o : opts, parse(",")
            gettoken comma o : o, parse(",")
        }
        use `base', clear
        generate long w = .
        * strata() leaves "" out, cluster() keeps it
        local empty = cond(strpos("`o'", "strata(g)"), "empty", "")
        _s2s_rankkey g rg `empty'
        _s2s_rankkey note rn
        local ro : subinstr local o "(g)" "(rg)", all
        local ro : subinstr local ro "(note)" "(rn)", all
        set seed 31
        cbsample `n', `ro'
        drop rg rn
        tempfile expected
        quietly save `expected'
        use `base', clear
        generate long w = .
        set seed 31
        cbsample `n', `o'
        cf _all using `expected'
    }
    if _rc test_fail "BSAMP-01 cbsample `opts'" "r(`=_rc')"
    else test_pass "BSAMP-01 cbsample `opts'"
}
* bsample's sizes: each stratum (strata() alone) keeps its count
capture noisily {
    use `base', clear
    generate long w = .
    set seed 3
    cbsample, strata(g) weight(w)
    _s2s_rankkey g rg
    tempfile got
    collapse (sum) w (count) id if g != "", by(rg)
    assert w == id
    use `base', clear
    generate long w = .
    set seed 3
    bsample, strata(g) weight(w)
    _s2s_rankkey g rg
    collapse (sum) w (count) id if g != "", by(rg)
    assert w == id
}
if _rc test_fail "BSAMP-01 cbsample strata(strL) sizes match bsample" "r(`=_rc')"
else test_pass "BSAMP-01 cbsample strata(strL) sizes match bsample"

* ---------------------------------------------------------------------------
* DESTR-06: cdestring with long and binary strL values is destring
* ---------------------------------------------------------------------------
capture program drop _s2s_destring_data
program define _s2s_destring_data
    clear
    quietly set obs 40
    tempfile b1 b2
    tempname fh
    file open `fh' using `b1', write binary replace
    file write `fh' %1bu (49) %1bu (50) %1bu (0)
    file close `fh'
    file open `fh' using `b2', write binary replace
    file write `fh' %1bu (0) %1bu (55)
    file close `fh'
    * a: numeric text, long padding and leading zeros
    generate strL a = string(_n * 3)
    quietly replace a = " " * 3000 + "5" in 3
    quietly replace a = "0" * 3000 + "12" in 4
    quietly replace a = "7" + " " * 9000 in 5
    quietly replace a = "." in 6
    quietly replace a = ".b" + " " * 2500 in 7
    quietly replace a = "" in 8
    * b: binary values, read as Stata reads them
    generate strL b = string(_n)
    quietly replace b = fileread("`b1'") in 2
    quietly replace b = fileread("`b2'") in 9
    * c: nonnumeric and too large
    generate strL c = string(_n / 8)
    quietly replace c = "1" * 400 in 10
    quietly replace c = "1" * 300 in 11
    quietly replace c = "x" + "1" * 3000 in 12
    * d: percent, commas and characters to ignore in long text
    generate strL d = string(_n) + "%"
    quietly replace d = "1,5%" + " " * 3000 in 13
    quietly replace d = "$" + "1" * 2 + "," * 2500 + "%" in 14
    * invalid UTF-8 after long text
    quietly replace d = "A" * 3000 + "4" + char(200) in 15
    * e: a fixed string, f: a strL the plugin reads
    generate str8 e = string(_n * 11)
    generate strL f = string(_n) + ".25"
    label variable a "label a"
    char a[mine] kept
end
_s2s_destring_data
tempfile dbase
quietly save `dbase'
local vars a b c d e f
foreach opts in "" "force" "float force" "percent force" "dpcomma force" ///
    `"ignore(" ") force"' `"ignore(",%$") force"' `"ignore(",%$") percent float force"' ///
    `"ignore("0") force"' `"ignore("A", illegal) force"' `"ignore("1 ", asbytes) force"' {
    capture noisily {
        use `dbase', clear
        destring `vars', generate(n1 n2 n3 n4 n5 n6) `opts'
        tempfile expected
        quietly save `expected'
        use `dbase', clear
        cdestring `vars', generate(n1 n2 n3 n4 n5 n6) `opts'
        cf _all using `expected'
        * the same variables, created or not, in the same places and types
        local got
        foreach v of varlist _all {
            local got `got' `v':`: type `v''
        }
        use `expected', clear
        local want
        foreach v of varlist _all {
            local want `want' `v':`: type `v''
        }
        assert "`got'" == "`want'"
        * replace
        use `dbase', clear
        destring `vars', replace `opts'
        quietly save `expected', replace
        use `dbase', clear
        cdestring `vars', replace `opts'
        cf _all using `expected'
        foreach v of local vars {
            local t1 : type `v'
            preserve
            use `expected', clear
            local t2 : type `v'
            restore
            assert "`t1'" == "`t2'"
        }
    }
    if _rc test_fail `"DESTR-06 cdestring long/binary strL `opts'"' "r(`=_rc')"
    else test_pass `"DESTR-06 cdestring long/binary strL `opts'"'
}
* if/in: the selected observations, others missing
capture noisily {
    use `dbase', clear
    keep in 1/8
    drop c d
    destring a b e f, generate(n1 n2 n3 n4)
    tempfile expected
    quietly save `expected'
    use `dbase', clear
    keep in 1/8
    drop c d
    cdestring a b e f if _n != 4 in 2/8, generate(n1 n2 n3 n4)
    foreach v in n1 n2 n3 n4 {
        rename `v' got_`v'
    }
    quietly merge 1:1 _n using `expected', nogenerate
    foreach v in n1 n2 n3 n4 {
        assert got_`v' == cond(_n != 4 & _n >= 2, `v', .)
    }
}
if _rc test_fail "DESTR-06 cdestring long/binary strL with if/in" "r(`=_rc')"
else test_pass "DESTR-06 cdestring long/binary strL with if/in"

* ---------------------------------------------------------------------------
* RJOIN-01: crangejoin with strL master and using variables is rangejoin
* ---------------------------------------------------------------------------
use `base', clear
keep id x note h s
rename note unote
generate strL ulong = "u" * 2500 + string(_n)
rename id uid
label variable ulong "using long"
note ulong: using note
tempfile ru
quietly save `ru'
foreach spec in "x -2 2 using `ru'" "x -2 2 using `ru', by(s)" ///
    "x -1 1 using `ru', keepusing(unote uid)" "x 0 0 using `ru', all prefix(u_)" ///
    "x 50 60 using `ru'" "x . 3 using `ru', by(s) keepusing(ulong h)" {
    capture noisily {
        use `base', clear
        drop g
        rangejoin `spec'
        tempfile expected
        quietly save `expected'
        use `base', clear
        drop g
        crangejoin `spec'
        cf _all using `expected'
        unab got : _all
        use `expected', clear
        unab want : _all
        assert "`got'" == "`want'"
        foreach v of local want {
            local t1 : type `v'
            local f1 : format `v'
            local l1 : variable label `v'
            preserve
            use `expected', clear
            assert "`t1'" == "`: type `v''" & "`f1'" == "`: format `v''"
            assert `"`l1'"' == `"`: variable label `v''"'
            restore
        }
    }
    if _rc test_fail "RJOIN-01 crangejoin `spec'" "r(`=_rc')"
    else test_pass "RJOIN-01 crangejoin `spec'"
}
* snapshot and fast give the same join
foreach mode in snapshot fast {
    capture noisily {
        use `base', clear
        drop g
        rangejoin x -2 2 using `ru', by(s)
        tempfile expected
        quietly save `expected'
        use `base', clear
        drop g
        crangejoin x -2 2 using `ru', by(s) `mode'
        cf _all using `expected'
    }
    if _rc test_fail "RJOIN-01 crangejoin strL variables, `mode'" "r(`=_rc')"
    else test_pass "RJOIN-01 crangejoin strL variables, `mode'"
}
* notes and characteristics of master strL variables survive
capture noisily {
    use `base', clear
    drop g
    crangejoin x -2 2 using `ru', by(s)
    assert `"`: char note[note1]'"' == "a note on the payload"
    assert `"`: variable label note'"' == "payload"
    assert `"`: char ulong[note1]'"' == "using note"
}
if _rc test_fail "RJOIN-01 crangejoin keeps strL labels and notes" "r(`=_rc')"
else test_pass "RJOIN-01 crangejoin keeps strL labels and notes"
* strL by() keys: rangejoin's merge refuses them, with merge's messages
foreach side in master using both {
    capture noisily {
        use `ru', clear
        generate strL key2 = "k" + string(mod(_n, 3))
        if "`side'" == "master" quietly recast str2 key2
        tempfile ru2
        quietly save `ru2'
        use `base', clear
        generate strL key2 = "k" + string(mod(_n, 3))
        if "`side'" == "using" quietly recast str2 key2
        tempfile m2
        quietly save `m2'
        tempfile log1 log2
        quietly log using `log1', text replace name(s2s1)
        capture noisily rangejoin x -1 1 using `ru2', by(key2)
        local rc1 = _rc
        quietly log close s2s1
        cf _all using `m2'
        quietly log using `log2', text replace name(s2s2)
        capture noisily crangejoin x -1 1 using `ru2', by(key2)
        local rc2 = _rc
        quietly log close s2s2
        cf _all using `m2'
        assert `rc1' == 106 & `rc2' == 106
        * the same error lines (rangejoin also names rangestat's version)
        mata: _s2s_loglines(st_local("log1"), st_local("log2"))
    }
    if _rc test_fail "RJOIN-01 crangejoin by(strL) in `side' refused as rangejoin" "r(`=_rc')"
    else test_pass "RJOIN-01 crangejoin by(strL) in `side' refused as rangejoin"
}
* rollback with strL variables on both sides
use `base', clear
drop g
foreach fault in schema store post strl {
    local rc = cond("`fault'" == "store", 920, 1)
    _ctools_rb_case, cmd(crangejoin x -3 3 using `ru', by(s)) ///
        global(CTOOLS_crangejoin_fault) fault(`fault') rc(`rc') ///
        testname("RJOIN-01 rollback with strL variables: `fault'")
}

* ---------------------------------------------------------------------------
* MERGE-08: cmerge writes strL values (updates, using-only rows) as merge
* ---------------------------------------------------------------------------
use `base', clear
keep if mod(id, 3)
keep id s note y
quietly replace note = "" if mod(id, 5) == 0
quietly replace y = . if mod(id, 7) == 0
tempfile mm
quietly save `mm'
use `base', clear
keep if mod(id, 2)
keep id s note g y
quietly replace id = id + 40
quietly replace y = y * 10
tempfile mu
quietly save `mu'
foreach spec in "1:1 id using `mu'" "1:1 id using `mu', update" ///
    "1:1 id using `mu', update replace" "1:1 id using `mu', keepusing(g)" ///
    "1:1 id using `mu', keep(master match) update" "m:1 s using `mu'b" ///
    "1:1 id using `mu', nogenerate keepusing(note) update replace" {
    capture noisily {
        if strpos("`spec'", "`mu'b") {
            use `mu', clear
            bysort s (id): keep if _n == 1
            tempfile mub
            quietly save `mub'
            local spec : subinstr local spec "`mu'b" "`mub'"
        }
        * merge breaks ties in its sort at random: the same seed for both
        use `mm', clear
        set sortseed 4021
        merge `spec'
        tempfile expected
        quietly save `expected'
        use `mm', clear
        set sortseed 4021
        cmerge `spec'
        cf _all using `expected'
        unab got : _all
        foreach v of local got {
            local t1 : type `v'
            preserve
            use `expected', clear
            assert "`t1'" == "`: type `v''"
            restore
        }
    }
    if _rc test_fail "MERGE-08 cmerge `spec'" "r(`=_rc')"
    else test_pass "MERGE-08 cmerge `spec'"
}
* a strL key: merge's error
capture noisily {
    use `mu', clear
    generate strL key = string(id)
    tempfile mukey
    quietly save `mukey'
    use `mm', clear
    generate strL key = string(id)
    capture merge 1:1 key using `mukey'
    local rc1 = _rc
    capture cmerge 1:1 key using `mukey'
    assert _rc == `rc1' & `rc1' == 106
}
if _rc test_fail "MERGE-08 cmerge strL key refused as merge" "r(`=_rc')"
else test_pass "MERGE-08 cmerge strL key refused as merge"

* ---------------------------------------------------------------------------
* MERGE-09: cmerge, tempdta writes merges with strL variables through a .dta
* and use, as merge leaves them: values, storage types, formats, labels,
* value labels, notes, characteristics, data label, sort marker, c(changed)
* ---------------------------------------------------------------------------
capture mata: mata drop _cm9_meta()
capture mata: mata drop _cm9_keep()
capture mata: mata drop _cm9_same()
mata:
string scalar _cm9_meta(string scalar extra)
{
    string scalar out, owner
    string colvector cn, labs, text
    real colvector values
    real scalar j, i

    out = extra
    for (j = 1; j <= st_nvar(); j++) {
        out = out + "|" + st_varname(j) + ":" + st_vartype(j) + ":" + st_varformat(j) + ":" +
              st_varlabel(j) + ":" + st_varvaluelabel(j)
    }
    for (j = 0; j <= st_nvar(); j++) {
        owner = j ? st_varname(j) : "_dta"
        cn = st_dir("char", owner, "*")
        if (rows(cn)) cn = sort(cn, 1)
        for (i = 1; i <= rows(cn); i++) {
            out = out + "|" + owner + "[" + cn[i] + "]=" + st_global(owner + "[" + cn[i] + "]")
        }
    }
    stata("quietly label dir")
    labs = tokens(st_global("r(names)"))'
    if (rows(labs)) labs = sort(labs, 1)
    for (i = 1; i <= rows(labs); i++) {
        st_vlload(labs[i], values, text)
        out = out + "|L:" + labs[i] + ":" + invtokens(strofreal(values', "%21.0g")) + ":" + invtokens(text')
    }
    return(out)
}

/* The signature in external __cm9_<name>; compared in Mata, since note
   texts and labels can hold quotes that a macro cannot carry */
void _cm9_keep(string scalar name, string scalar extra)
{
    pointer() scalar p

    rmexternal("__cm9_" + name)
    p = crexternal("__cm9_" + name)
    *p = _cm9_meta(extra)
}

void _cm9_same(string scalar a, string scalar b)
{
    if (*findexternal("__cm9_" + a) != *findexternal("__cm9_" + b)) exit(error(9))
}
end
capture program drop _cm9_sig
program define _cm9_sig
    args name
    local dl : data label
    local sb : sortedby
    mata: _cm9_keep("`name'", st_local("dl") + "#" + st_local("sb") + "#" + strofreal(c("changed")))
end

* master and using with labels, notes, characteristics and data labels; the
* using file defines bl differently (the master's definition stays) and an
* unattached label of its own. Both hold the same note on note (merge copies
* only the using notes the master does not have).
use `base', clear
keep if mod(id, 3)
keep id s note y
quietly replace note = "" if mod(id, 5) == 0
quietly replace y = . if mod(id, 7) == 0
generate byte b = 1 + mod(id, 2)
label define bl 1 "one" 2 "two"
label define unused 7 "seven"
label values b bl
label variable note "master note"
note note: master's note on note
note: master's dataset note
char _dta[mchar] "m"
char y[ych] "y characteristic"
label data "master data"
* in key order with ties broken by id, so merge's order within a key is
* fixed (merge breaks ties at random when it sorts)
sort s id
tempfile mm9
quietly save `mm9'
use `base', clear
keep if mod(id, 2)
keep id s note g y
quietly replace id = id + 40
quietly replace y = y * 10
generate byte z = 1 + mod(id, 3)
label define zl 1 "a" 2 "b" 3 "c"
label define bl 1 "uno" 2 "dos"
label values z zl
label variable g "using g"
note g: using's note on g
note note: using's note on note
char _dta[uchar] "u"
label data "using data"
sort s id
tempfile mu9 mub9 mmb9
quietly save `mu9'
bysort s (id): keep if _n == 1
quietly save `mub9'
use `mm9', clear
bysort s (id): keep if _n == 1
quietly save `mmb9'

local specs `""1:1 id using `mu9'" "1:1 id using `mu9', update""'
local specs `"`specs' "1:1 id using `mu9', update replace""'
local specs `"`specs' "1:1 id using `mu9', keepusing(g z)""'
local specs `"`specs' "1:1 id using `mu9', keep(master match) update""'
local specs `"`specs' "1:1 id using `mu9', nogenerate keepusing(note) update replace""'
local specs `"`specs' "1:1 id using `mu9', generate(how) nolabel nonotes""'
local specs `"`specs' "1:1 id using `mu9', keep(using)""'
local specs `"`specs' "m:1 s using `mub9'" "m:1 s using `mub9', keepusing(note) update""'
local specs `"`specs' "1:1 _n using `mu9'" "m:m s using `mu9'""'
local i 0
foreach spec of local specs {
    local ++i
    if strpos(`"`spec'"', ",") local with `"`spec' tempdta"'
    else local with `"`spec', tempdta"'
    local master `mm9'
    if strpos(`"`spec'"', "1:m") local master `mmb9'
    capture noisily {
        use `master', clear
        set sortseed 4021
        merge `spec'
        _cm9_sig want
        tempfile expected
        quietly save `expected'
        use `master', clear
        set sortseed 4021
        cmerge `with'
        * the plugin route, not merge itself (which sets no r(time))
        local rs : r(scalars)
        assert `: list posof "time" in rs'
        _cm9_sig got
        mata: _cm9_same("want", "got")
        cf _all using `expected'
    }
    if _rc test_fail "MERGE-09 cmerge spec `i' with tempdta" "r(`=_rc')"
    else test_pass "MERGE-09 cmerge spec `i' with tempdta"
}

* 1:m (a master unique in s), and a merge without strL variables (tempdta
* changes nothing)
foreach case in "1:m s using `mu9'" "1:1 id using `mu9', keepusing(z y)" {
    capture noisily {
        use `mmb9', clear
        if strpos(`"`case'"', "keepusing(z y)") {
            use `mm9', clear
            drop note
        }
        tempfile start
        quietly save `start'
        set sortseed 4021
        merge `case'
        _cm9_sig want
        tempfile expected
        quietly save `expected'
        use `start', clear
        set sortseed 4021
        if strpos(`"`case'"', ",") cmerge `case' tempdta
        else cmerge `case', tempdta
        _cm9_sig got
        mata: _cm9_same("want", "got")
        cf _all using `expected'
    }
    local how = cond(strpos(`"`case'"', "1:m"), "1:m", "without strL")
    if _rc test_fail "MERGE-09 cmerge `how' tempdta" "r(`=_rc')"
    else test_pass "MERGE-09 cmerge `how' tempdta"
}

* 30,000 master rows, m:1, many blocks of the writer
capture noisily {
    clear
    quietly set obs 30000
    generate long k = mod(_n * 7, 900)
    generate strL t = "row " + string(_n) + (100 * "x") * (mod(_n, 4) == 0)
    generate double x = _n / 3
    * in key order (merge breaks ties at random when it sorts)
    sort k x
    tempfile bigm bigu
    quietly save `bigm'
    clear
    quietly set obs 1000
    generate long k = _n - 1
    generate strL u = "key " + string(k) + (3000 * "y") * (mod(k, 50) == 0)
    quietly save `bigu'
    use `bigm', clear
    merge m:1 k using `bigu'
    tempfile expected
    quietly save `expected'
    use `bigm', clear
    cmerge m:1 k using `bigu', tempdta
    cf _all using `expected'
}
if _rc test_fail "MERGE-09 cmerge m:1 30,000 rows tempdta" "r(`=_rc')"
else test_pass "MERGE-09 cmerge m:1 30,000 rows tempdta"

* A failure or Break before use leaves the master as it was; after use the
* merged data stand (there is nothing left to undo)
foreach fault in schema store post dtaused {
    capture noisily {
        use `mm9', clear
        _cm9_sig before
        quietly datasignature
        local ds `"`r(datasignature)'"'
        global CTOOLS_cmerge_fault `fault'
        capture cmerge 1:1 id using `mu9', tempdta
        local rc = _rc
        global CTOOLS_cmerge_fault
        assert `rc' == cond("`fault'" == "store", 920, 1)
        if "`fault'" != "dtaused" {
            _cm9_sig after
            mata: _cm9_same("before", "after")
            quietly datasignature
            assert `"`r(datasignature)'"' == `"`ds'"'
        }
        else {
            use `mm9', clear
            merge 1:1 id using `mu9'
            tempfile expected
            quietly save `expected'
            use `mm9', clear
            global CTOOLS_cmerge_fault dtaused
            capture cmerge 1:1 id using `mu9', tempdta
            global CTOOLS_cmerge_fault
            cf _all using `expected'
        }
    }
    local rc = _rc
    global CTOOLS_cmerge_fault
    if `rc' test_fail "MERGE-09 cmerge tempdta fault `fault'" "r(`rc')"
    else test_pass "MERGE-09 cmerge tempdta fault `fault'"
}

print_summary s2_strl
if $TESTS_FAILED > 0 exit 9
global CTOOLS_COMPONENT_COMPLETE "s2_strl"
    global CTOOLS_CASE_COMPLETE "s2_strl"
    exit
}

di as error "Unknown s2_strl case: `__case'"
exit 198
