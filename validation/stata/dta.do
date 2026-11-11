* dta: the plugin's .dta writer (internal plugin command "dta", used by cimport).
* Named cases share this file, each do-call has its own local scope.
gettoken __case 0 : 0
macro shift
if "`__case'" == "" local __case "dta"

if "`__case'" == "dta" {
do "validation/support/validate_setup.do"
capture mkdir "temp"
_ctools_load
capture program ctools_plugin, plugin using("`__ctools_plugin'")

* Write the dataset in memory through the writer (dta save) to `using',
* carrying each variable's storage type, format, variable label and value
* label name; threads() and sort() pass through.
capture program drop _dta_save
program define _dta_save
    syntax using/ [, threads(integer 0) sort(string)]
    local __ctools_dta_path `"`using'"'
    local __ctools_dta_types
    local __ctools_dta_names
    local __ctools_dta_formats
    local __ctools_dta_vallabs
    local j 0
    foreach v of varlist _all {
        local ++j
        local __ctools_dta_types `__ctools_dta_types' `: type `v''
        local __ctools_dta_names `__ctools_dta_names' `v'
        local __ctools_dta_formats `__ctools_dta_formats' `: format `v''
        local vl : value label `v'
        if "`vl'" == "" local vl "."
        local __ctools_dta_vallabs `__ctools_dta_vallabs' `vl'
        local __ctools_dta_label`j' `"`: variable label `v''"'
    }
    local __ctools_dta_datalabel `"`: data label'"'
    local __ctools_dta_threads `threads'
    local __ctools_dta_sort
    foreach v of local sort {
        unab v : `v'
        local j 0
        foreach w of varlist _all {
            local ++j
            if "`w'" == "`v'" local __ctools_dta_sort `__ctools_dta_sort' `j'
        }
    }
    if c(k) plugin call ctools_plugin _all, "dta save"
    else plugin call ctools_plugin, "dta save"
end

* The same dataset with every attribute the save carries: storage types,
* formats, variable labels, value label names, data label and values.
capture program drop _dta_same
program define _dta_same
    syntax using/
    preserve
    quietly describe, replace clear
    tempfile mine
    quietly save `mine'
    restore
    preserve
    quietly use `"`using'"', clear
    quietly describe, replace clear
    cf _all using `mine'
    restore
    cf _all using `"`using'"'
end

* ----------------------------------------------------------------- selftest
* Every storage type, the 27 missing values, boundary values, str1 and a
* full-width str2045, strLs (empty, 1, 2,045, 2,046 and 100,000 bytes,
* binary with and without NUL, Unicode), a value label given out of order,
* notes, characteristics, a data label and a sort marker; formats 118 and
* 119, one, two and the default number of threads, and 1, 3,000 and 50,000
* rows (many blocks and waves).
capture program drop _dta_check_selftest
program define _dta_check_selftest
    args N strls
    assert _N == `N'
    local types "id long b byte i int l long f float d double s1 str1 s2045 str2045 L strL L2 strL"
    while "`types'" != "" {
        gettoken v types : types
        gettoken t types : types
        assert "`: type `v''" == "`t'"
    }
    assert "`: format d'" == "%12.4f"
    assert "`: format id'" == "%12.0g" & "`: format b'" == "%8.0g" & "`: format f'" == "%9.0g"
    assert "`: format s1'" == "%9s" & "`: format s2045'" == "%2045s" & "`: format L'" == "%9s"
    assert "`: data label'" == "ctools dta selftest"
    assert "`: sortedby'" == "id"
    assert "`: value label b'" == "yn"
    assert `"`: label yn 1'"' == "one" & `"`: label yn -1'"' == "minus one" & `"`: label yn 0'"' == "zero"
    assert `"`: variable label d'"' == "a double" & `"`: variable label L'"' == "a strL"
    assert `"`: char _dta[note1]'"' == "dataset note" & `"`: char d[note1]'"' == "double note"
    assert `"`: char s1[mychar]'"' == "char text"
    tempvar r m k
    quietly generate long `r' = _n - 1
    quietly generate byte `m' = mod(`r', 32)
    quietly generate byte `k' = `m' - 27
    assert id == `r' + 1
    forvalues j = 0/26 {
        local mv = cond(`j' == 0, ".", "." + char(96 + `j'))
        foreach v in b i l f d {
            assert `v' == `mv' if `m' == `j'
        }
    }
    assert b == cond(`k' == 0, -127, cond(`k' == 1, 100, cond(`k' == 2, 0, cond(`k' == 3, -1, 1)))) if `m' >= 27
    assert i == cond(`k' == 0, -32767, cond(`k' == 1, 32740, cond(`k' == 2, 0, cond(`k' == 3, -1, 1)))) if `m' >= 27
    assert l == cond(`k' == 0, -2147483647, cond(`k' == 1, 2147483620, cond(`k' == 2, 0, cond(`k' == 3, -1, 1)))) if `m' >= 27
    * c(maxfloat) is a rounded double; float() gives the largest float
    assert f == cond(`k' == 0, -float(c(maxfloat)), cond(`k' == 1, float(c(maxfloat)), cond(`k' == 2, 0.5, cond(`k' == 3, -0.25, 3)))) if `m' >= 27
    assert d == cond(`k' == 0, c(mindouble), cond(`k' == 1, c(maxdouble), cond(`k' == 2, 0.1, cond(`k' == 3, -1e-300, 1)))) if `m' >= 27
    assert s1 == cond(mod(`r', 3) == 0, "", cond(mod(`r', 3) == 1, "a", "z"))
    assert s2045 == char(65 + mod(`r', 26)) * 2045 if mod(`r', 2)
    assert s2045 == "short" + string(`r') if !mod(`r', 2)
    quietly replace `m' = mod(`r', 7)
    assert L == "" if `m' == 0
    assert L == "x" if `m' == 1
    assert L == 2045 * "b" if `m' == 2
    assert L == 2046 * "c" if `m' == 3
    assert L == 100000 * "d" if `m' == 4
    assert L == "a" + char(0) + "b" & strlen(L) == 3 if `m' == 5
    assert L == "héllo wörld" if `m' == 6
    assert L2 == cond(mod(`r', 5), "", "five" + string(`r'))
    quietly count if L != ""
    local n = r(N)
    quietly count if L2 != ""
    assert `strls' == `n' + r(N)
end

tempfile st
local cases ""
foreach rel in 118 119 {
    foreach t in 0 1 2 {
        foreach N in 1 3000 {
            local cases `cases' `rel':`t':`N'
        }
    }
}
local cases `cases' 118:0:50000 119:2:50000 118:1:50000
foreach c of local cases {
    tokenize "`c'", parse(":")
    local rel `1'
    local t `3'
    local N `5'
    capture noisily {
        clear
        local __ctools_dta_path "`st'"
        plugin call ctools_plugin, "dta selftest `rel' `t' `N'"
        local strls `__ctools_dta_strls'
        * the release Stata reads in the header
        mata: fh = fopen(st_local("st"), "r"); st_local("head", fread(fh, 31)); fclose(fh)
        assert substr("`head'", 29, 3) == "`rel'"
        use "`st'", clear
        _dta_check_selftest `N' `strls'
        * Stata's own save of what was read gives back the same data
        tempfile again
        quietly save `again'
        _dta_same using `again'
    }
    if _rc test_fail "DTA-01 selftest release `rel' threads `t' N `N'" "r(`=_rc')"
    else test_pass "DTA-01 selftest release `rel' threads `t' N `N'"
}

* ------------------------------------------------------------- round trip
* A random dataset through dta save, compared with the original by cf and
* by its description: types, formats, labels and value label names. Callers
* set their own seeds; this one only seeds the case before the program that
* draws is defined, as the seed check requires.
set seed 20261109
capture program drop _dta_random
program define _dta_random
    args N
    clear
    quietly set obs `N'
    tempfile b1 b2
    tempname fh
    file open `fh' using `b1', write binary replace
    file write `fh' %1bu (0) %1bu (65) %1bu (0) %1bu (66)
    file close `fh'
    file open `fh' using `b2', write binary replace
    file write `fh' %1bu (200) %1bu (49) %1bu (255) %1bu (32)
    file close `fh'
    generate byte b = floor(runiform() * 228) - 127
    generate int i = floor(runiform() * 65508) - 32767
    generate long l = floor(runiform() * 4294967268) - 2147483647
    generate float f = rnormal() * 1e10
    generate double d = rnormal() * 10^floor(runiform() * 600 - 300)
    foreach v in b i l f d {
        quietly replace `v' = . if runiform() < .05
        quietly replace `v' = .a if runiform() < .02
        quietly replace `v' = .z if runiform() < .02
        quietly replace `v' = .m if runiform() < .02
    }
    generate str7 s = substr("abcdefghij", 1, floor(runiform() * 8))
    generate str244 w = 244 * "w" if runiform() < .5
    generate strL L = ""
    quietly replace L = "x" * floor(runiform() * 5000) if runiform() < .6
    quietly replace L = fileread("`b1'") if runiform() < .05
    quietly replace L = fileread("`b2'") if runiform() < .05
    quietly replace L = "Ωmega " + string(_n) if runiform() < .1
    label define lab 1 "one" 2 "two"
    label values b lab
    label variable d "a double, labeled"
    label variable L "strL"
    format d %20.6e
    format s %-10s
    label data "random data"
end

foreach N in 0 1 997 40000 {
    foreach t in 0 1 2 {
        capture noisily {
            set seed `=1000 + `N' + `t''
            _dta_random `N'
            tempfile orig mine
            quietly save `orig'
            _dta_save using `mine', threads(`t')
            use `mine', clear
            assert "`: data label'" == "random data"
            assert "`: sortedby'" == ""
            _dta_same using `orig'
        }
        if _rc test_fail "DTA-02 round trip N `N' threads `t'" "r(`=_rc')"
        else test_pass "DTA-02 round trip N `N' threads `t'"
    }
}

* a sort marker set by the caller, and a dataset of observations only
capture noisily {
    set seed 77
    _dta_random 500
    sort i b
    tempfile orig mine
    quietly save `orig'
    _dta_save using `mine', sort(i b)
    use `mine', clear
    assert "`: sortedby'" == "i b"
    _dta_same using `orig'
    clear
    quietly set obs 5
    _dta_save using `mine'
    use `mine', clear
    assert _N == 5 & c(k) == 0
}
if _rc test_fail "DTA-03 sort marker and no variables" "r(`=_rc')"
else test_pass "DTA-03 sort marker and no variables"

* ---------------------------------------------------------------- refusals
* Specs the writer refuses, leaving no file behind
local expect "dupname:198 badname:198 reserved:198 strformat:198 byte:459 fraction:459 strwidth:459"
foreach c of local expect {
    tokenize "`c'", parse(":")
    capture noisily {
        tempfile bad
        local __ctools_dta_path "`bad'"
        plugin call ctools_plugin, "dta error `1'"
        assert "`__ctools_dta_rc'" == "`3'"
        assert `"`__ctools_dta_message'"' != ""
        capture confirm file "`bad'"
        assert _rc == 601
    }
    if _rc test_fail "DTA-05 refuses `1' with r(`3')" "r(`=_rc') got `__ctools_dta_rc': `__ctools_dta_message'"
    else test_pass "DTA-05 refuses `1' with r(`3')"
}

* ---------------------------------------------------------------- producer
* Rows from a producer (groups of 1,000) next to a column in memory. keep 1
* keeps no group, so each is filled again for its strLs; 0 keeps them all.
capture program drop _dta_check_producer
program define _dta_check_producer
    args N
    assert _N == `N'
    assert "`: type pid'" == "long" & "`: type pnum'" == "double" & "`: type pstr'" == "str8"
    assert "`: type pL'" == "strL" & "`: type m'" == "double"
    tempvar r
    quietly generate long `r' = _n - 1
    assert pid == `r' + 1
    assert pnum == .b if !mod(`r', 13)
    assert pnum == `r' * 0.25 if mod(`r', 13)
    assert pstr == "p" + string(mod(`r', 1000))
    assert pL == "" if mod(`r', 4) == 0
    assert pL == char(97 + mod(`r', 26)) * 3000 if mod(`r', 4) == 1
    assert pL == "short" + string(`r') if mod(`r', 4) == 2
    assert pL == char(0) + string(`r') if mod(`r', 4) == 3
    assert m == `r' * 2
end
foreach c in 2500:0:0 2500:1:1 2500:2:1 2500:0:1 20001:0:0 20001:2:1 {
    tokenize "`c'", parse(":")
    local N `1'
    local t `3'
    local keep `5'
    capture noisily {
        clear
        tempfile pf
        local __ctools_dta_path "`pf'"
        plugin call ctools_plugin, "dta producer `N' `t' `keep'"
        local groups = ceil(`N' / 1000)
        assert `__ctools_dta_fills' == `__ctools_dta_releases'
        assert `__ctools_dta_fills' == cond(`keep' == 1, 2, 1) * `groups'
        use `pf', clear
        _dta_check_producer `N'
    }
    if _rc test_fail "DTA-06 producer N `N' threads `t' keep `keep'" "r(`=_rc')"
    else test_pass "DTA-06 producer N `N' threads `t' keep `keep'"
}

* ------------------------------------------------------------ cimport strL
* With option tempdta, an import that may replace the data and whose strL
* text outweighs its other columns arrives through the writer: the whole
* dataset (cimport fulldta), or, when that declines, its strL columns
* (cimport strldta). Without tempdta, the strL columns travel through Mata. Both must match import delimited. Quoted
* values with doubled quotes and commas, empty and Unicode values, two strL
* columns and 70,000 rows (many groups); the short values are long enough
* that the strL text outweighs id and num.
capture noisily {
    tempfile csv native
    tempname fh
    file open `fh' using "`csv'", write text replace
    file write `fh' "id,long1,num,long2" _n
    local pad = 40 * "p"
    forvalues r = 1/70000 {
        if `r' == 5 local a = `"""' + 1200 * "q" + `"""",""' + 1200 * "z" + `"""'
        else if `r' == 7 local a = 3000 * "a"
        else if `r' == 9 local a `""line1, Ωé ""quoted""""'
        else if `r' == 69999 local a = 2100 * "e"
        else if mod(`r', 3) == 0 local a ""
        else local a "v`r'`pad'"
        local b = cond(mod(`r', 1000) == 1, 2046 * "b", cond(mod(`r', 2), "", "w`r'`pad'"))
        file write `fh' `"`r',`a',`=`r' / 4',`b'"' _n
    }
    file close `fh'
    clear
    import delimited using "`csv'", clear encoding(utf-8)
    quietly save `native'
}
local made = _rc
* default: Mata. fulldta: option tempdta. strldta: tempdta with the ado's
* validation hook making fulldta decline, so the strL columns come as a .dta
* and the other variables are added after.
foreach how in default fulldta strldta {
    if "`how'" == "strldta" global CTOOLS_cimport_fault "fulldta"
    capture noisily {
        assert `made' == 0
        clear
        cimport delimited using "`csv'", clear encoding(utf-8) `=cond("`how'" == "default", "", "tempdta")'
        assert "`: type long1'" == "strL" & "`: type long2'" == "strL"
        assert c(changed) == 1
        cf _all using `native'
    }
    local rc = _rc
    global CTOOLS_cimport_fault
    if `rc' test_fail "DTA-07 cimport strL columns, `how'" "r(`rc')"
    else test_pass "DTA-07 cimport strL columns, `how'"
}

* ------------------------------------------------------------------ wide
* More than 32,767 variables: format 119 (last: it clears Mata)
capture noisily {
    if c(MP) & c(stata_version) >= 15 {
        clear
        * set maxvar needs an empty Mata; nothing after this test uses Mata
        clear mata
        set maxvar 33000
        tempfile wide
        local __ctools_dta_path "`wide'"
        plugin call ctools_plugin, "dta wide 32768"
        mata: fh = fopen(st_local("wide"), "r"); st_local("head", fread(fh, 31)); fclose(fh)
        assert substr("`head'", 29, 3) == "119"
        use `wide', clear
        assert c(k) == 32768 & _N == 3
        assert v1 == _n - 1 & v32768 == mod(32767 + _n - 1, 100)
        assert "`: type v32768'" == "byte"
        clear
        clear mata
        set maxvar 5000
    }
}
if _rc test_fail "DTA-04 32,768 variables as format 119" "r(`=_rc')"
else test_pass "DTA-04 32,768 variables as format 119"

print_summary dta
if $TESTS_FAILED > 0 exit 9
global CTOOLS_COMPONENT_COMPLETE "dta"
    global CTOOLS_CASE_COMPLETE "dta"
    exit
}

di as error "Unknown dta case: `__case'"
exit 198
