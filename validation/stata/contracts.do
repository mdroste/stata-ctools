* contracts: named cases share this file, each do-call has its own local scope.
gettoken __case 0 : 0
macro shift
if "`__case'" == "" local __case "contracts"

if "`__case'" == "contracts" {
/* Boundary regressions: invalid wire arguments must not mutate observations. */
do validation/support/validate_setup.do
clear
set obs 3
gen double result = 17
_ctools_load
program p1_contracts_plugin, plugin using("`__ctools_plugin'")
foreach args in "1 0 count=18446744073709551616" "1 0 percent=nan" "2147483648 0 count=1" "1 -1 count=1" "1 0 count=1 count=2" "1 0 seedhi=1 count=1" {
    capture plugin call p1_contracts_plugin result, "csample `args'"
    local rc = _rc
    if `rc' == 198 test_pass "rejected malformed sampling arguments"
    else test_fail "malformed sampling arguments" "expected 198, got `rc'"
    assert result == 17
}
foreach cmd in "csplit write" "crangejoin prepare 1" "crangejoin write" "cmerge execute" "cio column 1" {
    capture plugin call p1_contracts_plugin result, "`cmd'"
    local rc = _rc
    if `rc' == 198 test_pass "rejected stale command phase"
    else test_fail "stale command phase" "expected 198, got `rc'"
    assert result == 17
}
* Failure recovery and repeated cleanup remain safe in the real dispatcher.
foreach cmd in "cimport clear" "cio clear" "cmerge clear" "csplit clear" "crangejoin clear" {
    plugin call p1_contracts_plugin, "`cmd'"
    plugin call p1_contracts_plugin, "`cmd'"
    test_pass "repeated cleanup"
}
csample 100
assert _N == 3 & result == 17
test_pass "successful command after malformed calls"
program drop p1_contracts_plugin
print_summary "contracts"
global CTOOLS_COMPONENT_COMPLETE "contracts"
    global CTOOLS_CASE_COMPLETE "contracts"
    exit
}

if "`__case'" == "bigdata" {
* args: directory containing the plugin and ado files, output log path
args plugin_dir
if `"`plugin_dir'"' == "" local plugin_dir "build"
clear all
set more off
adopath ++ "`plugin_dir'"
* Duplicates that do not match the opposite side; empty-side bypasses.
foreach side in master using {
    foreach empty in 0 1 {
        foreach strings in 0 1 {
            clear
            set obs 3
            gen long k = cond(_n == 1,1,2)
            gen str12 payload = "duplicate"
            if `strings' tostring k, replace
            tempfile duplicate other
            save `duplicate'
            clear
            set obs 1
            gen long k = 99
            gen str12 other = "unique"
            if `strings' tostring k, replace
            if `empty' keep if 0
            save `other'
            if "`side'" == "master" {
                use `duplicate', clear
                capture cmerge 1:1 k using `other'
            }
            else {
                use `other', clear
                capture cmerge 1:1 k using `duplicate'
            }
            assert _rc == 459
        }
    }
}
* Borrowed key strings, using-only keys, duplicates and shared payload.
clear
input str4 k double x str10 s
"b" 2 "bee"
"" 0 "empty"
"a" 1 "aye"
"b" 3 "bee2"
end
tempfile m u ref
save `m'
clear
input str4 k double z str10 t
"d" 4 "dee"
"a" 1 "aye"
"b" 2 "bee"
"" 0 "empty"
end
save `u'
use `m', clear
merge m:1 k using `u'
sort k x
save `ref'
use `m', clear
cmerge m:1 k using `u'
sort k x
cf _all using `ref'
* Sorted inputs with inserted using-only rows, shared columns, and duplicate
* master groups. Compare every value against merge; include both key types.
foreach strings in 0 1 {
    foreach sorted_master in 0 1 {
        foreach shared in 0 1 {
            clear
            set obs 3
            gen long k = cond(_n == 2, 1, 3)
            gen double x = cond(_n == 1, 30, cond(_n == 2, 10, .))
            gen str12 s = cond(_n == 1, "three", cond(_n == 2, "one", ""))
            if `strings' tostring k, replace
            if `sorted_master' sort k
            tempfile master using reference
            save `master'
            clear
            set obs 2
            gen long k = _n + 1
            gen double x = cond(_n == 1, 20, 33)
            gen str12 s = cond(_n == 1, "two", "using3")
            gen double z = 100 * k
            if !`shared' drop x s
            if `strings' tostring k, replace
            sort k
            save `using'
            foreach opts in "" "update" "update replace" {
                use `master', clear
                local genopt ""
                if "`opts'" != "" local genopt "nogenerate"
                merge m:1 k using `using', `opts' `genopt'
                sort k x s
                save `reference', replace
                use `master', clear
                cmerge m:1 k using `using', `opts' `genopt'
                sort k x s
                cf _all using `reference'
            }
        }
    }
}
di "BIGDATA_EDGE_PASS"
    global CTOOLS_CASE_COMPLETE "bigdata"
    exit
}

if "`__case'" == "command_optimizations" {
* Focused reference coverage for optimized kernels and wrapper fast paths.
do validation/support/validate_setup.do
quietly {
clear
set seed 43211
set obs 70001
gen long id=_n
gen double x=floor(runiform()*100)
gen double y=rnormal()
replace y=. if runiform()<.6
gen double g=floor(runiform()*100)
replace g=.a if mod(_n,97)==0
replace g=.z if mod(_n,101)==0
foreach groups in "g" "x" "y" "x y" {
    sort id
    ipolate y x, gen(ref) by(`groups') epolate
    sort id
    foreach threads in 1 4 12 {
        cipolate y x, gen(got) by(`groups') epolate threads(`threads')
        assert_var_equal got ref $DEFAULT_SIGFIGS "radix/alias by(`groups') threads(`threads')"
        assert id==_n
        drop got
    }
    drop ref
}
clear
set obs 40001
gen str2045 s=cond(mod(_n,3)==0, 150*"a::東京: "+" end  ", "  :aa::|--z;::end  ")
foreach trim in "" "notrim" {
    foreach limit in 1 5 200 {
        split s if mod(_n,7), gen(ref) parse("::" ":" "|" "--" ";" "東京") `trim' limit(`limit')
        local k=r(k_new)
        foreach threads in 1 4 12 {
            csplit s if mod(_n,7), gen(got) parse("::" ":" "|" "--" ";" "東京") `trim' limit(`limit') threads(`threads')
            assert r(k_new)==`k'
            forvalues j=1/`k' {
                assert got`j'==ref`j'
            }
            drop got*
            test_pass "packed tokens `trim' limit(`limit') threads(`threads')"
        }
        drop ref*
    }
}
* One master row expands across many output tiles; every row must retain ties.
clear
set obs 50001
gen long uid=_n
gen double key=mod(_n,31)
gen str16 text="using"+string(_n)
tempfile u m expected
save `u'
clear
set obs 1
gen double key=15
gen str12 master="one master"
save `m'
rangejoin key -20 20 using `u'
save `expected'
foreach threads in 1 4 12 {
    use `m', clear
    crangejoin key -20 20 using `u', threads(`threads')
    cf _all using `expected'
    test_pass "one-master tiled join threads(`threads')"
}
* Repeated restored master datasets retain metadata and stable using ties.
clear
set obs 40001
gen long mid=_n
gen double g=mod(_n,7)
replace g=.a if mod(_n,17)==0
replace g=.z if mod(_n,23)==0
gen double key=mod(_n,101)
label variable mid "Original master row"
format key %12.3f
save `u', replace
replace key=key+.5
save `m', replace
rangejoin key -1 1 using `u', by(g)
save `expected', replace
foreach threads in 1 4 12 {
    use `m', clear
    crangejoin key -1 1 using `u', by(g) threads(`threads')
    cf _all using `expected'
    assert "`: variable label mid'"=="Original master row"
    assert "`: format key'"=="%12.3f"
    test_pass "indexed extended-missing groups threads(`threads')"
}
}
print_summary "command_optimizations"
global CTOOLS_COMPONENT_COMPLETE "command_optimizations"
    global CTOOLS_CASE_COMPLETE "command_optimizations"
    exit
}

if "`__case'" == "data_io" {
* Minimal test for data I/O refactoring
clear all
adopath ++ "build"
set obs 10

* Test 1: cencode with 1 variable (string)
gen str10 city = ""
replace city = "NYC" in 1
replace city = "LA" in 2
replace city = "NYC" in 3
replace city = "Chicago" in 4
replace city = "LA" in 5
replace city = "Boston" in 6
replace city = "NYC" in 7
replace city = "Chicago" in 8
replace city = "LA" in 9
replace city = "Boston" in 10

di "About to call cencode..."
cencode city, gen(city_num) verbose
di "cencode done: " _rc
list city city_num in 1/5

* Test 2: cdecode with 1 variable (numeric)
di "About to call cdecode..."
cdecode city_num, gen(city_str) verbose
di "cdecode done: " _rc
list city_num city_str in 1/5

* Test 3: cdestring with 1 variable
clear
set obs 5
gen str10 x = ""
replace x = "1.5" in 1
replace x = "2.3" in 2
replace x = "3.7" in 3
replace x = "4.1" in 4
replace x = "5.9" in 5

di "About to call cdestring with generate..."
cdestring x, gen(x_num) verbose
di "cdestring generate done: " _rc
list x x_num

di "ALL TESTS PASSED"
    global CTOOLS_CASE_COMPLETE "data_io"
    exit
}

if "`__case'" == "memory_safety" {
/* Memory-audit integration contracts; native fault injection lives in
   test_memory_safety_native.py. Run through the machine's stata wrapper. */
do validation/support/validate_setup.do

capture noisily {
    clear
    set obs 120
    gen double x = _n + sin(_n)
    gen double y = 3 + 2*x + cos(_n)
    gen byte g = mod(_n, 3)
    gen double large_g = cond(g==0, -2^63, cond(g==1, 2^62, 1e30))
    quietly creghdfe y x, absorb(g)
    scalar slope = _b[x]
    quietly creghdfe y x, absorb(large_g)
    assert reldif(_b[x], scalar(slope)) < 1e-10
    replace large_g = g + .25
    quietly creghdfe y x, absorb(large_g)
    assert reldif(_b[x], scalar(slope)) < 1e-10

    gen double w = 300000000
    foreach vce in "" "vce(robust)" "vce(cluster g)" {
        quietly creghdfe y x [fw=w], absorb(g) `vce'
        assert e(N)==36000000000
        assert reldif(_b[x], scalar(slope)) < 1e-10
        assert _se[x]>0 & _se[x]<.
    }
}
if _rc test_fail "large FE labels and frequency weights" "r(`=_rc')"
else test_pass "large FE labels and frequency weights"

capture noisily {
    clear
    set obs 512
    gen long key = _n
    gen byte g = floor((_n-1)/64)
    gen double x = mod(_n, 10)
    crangestat (min) lo=x (max) hi=x (median) med=x (mean) avg=x, ///
        interval(key 0 0) by(g) threads(8)
    assert lo==x & hi==x & med==x & avg==x
    drop lo hi med avg
    crangestat (median) med=x (mean) avg=x, interval(key 0 0) threads(8)
    assert med==x & avg==x
    capture crangestat (mean) invalid=x, interval(key 0 0) threads(2147483647)
    assert _rc==198
}
if _rc test_fail "range scratch paths and thread bound" "r(`=_rc')"
else test_pass "range scratch paths and thread bound"

capture noisily {
    tempfile original selected
    save `original'
    set seed 142
    csample, count(20) by(g) threads(1)
    assert _N==160
    sort key
    save `selected'
    use `original', clear
    set seed 142
    csample, count(20) by(g) threads(8)
    sort key
    cf _all using `selected'
    use `original', clear
    csample if key<=100, count(20) threads(8)
    count if key<=100
    assert r(N)==20
    use `original', clear
    gen double freq=0
    cbsample 20 if key<=100, weight(freq) threads(8)
    quietly summarize freq, meanonly
    assert r(sum)==20
}
if _rc test_fail "sampling grouping and filtered maps" "r(`=_rc')"
else test_pass "sampling grouping and filtered maps"

capture noisily {
    clear
    set obs 1000
    set seed 324
    gen double x = rnormal()
    gen double z = x + rnormal()
    gen double z2 = 2*z
    gen double y = 2*x + 3*z + rnormal()
    gen g = mod(_n,10)
    foreach adjustment in "" "absorb(g)" {
        cbinscatter y x, controls(z) method(binsreg) nograph `adjustment'
        matrix reference=e(bindata)
        cbinscatter y x, controls(z z2) method(binsreg) nograph `adjustment'
        mata: assert(max(abs(st_matrix("reference")-st_matrix("e(bindata)"))) < 1e-7)
    }
}
if _rc test_fail "cbinscatter controls and FE offsets" "r(`=_rc')"
else test_pass "cbinscatter controls and FE offsets"

print_summary "memory safety"
global CTOOLS_COMPONENT_COMPLETE "memory_safety"
    global CTOOLS_CASE_COMPLETE "memory_safety"
    exit
}

if "`__case'" == "performance_paths" {
* Sample selection and complete destination writes for the October 2026 fast paths.
version 16
clear all
set more off
adopath ++ build
local candidate : environment CTOOLS_TEST_PLUGIN_DIR
if `"`candidate'"' != "" adopath ++ `"`candidate'"'
set seed 72391
set obs 20000
gen long id = _n
gen str8 g = "g" + string(mod(id,13))
gen double x = floor(id/13)
gen double y = 3*x + mod(id,13)
replace y = . if mod(id,5)==0
replace x = .a if mod(id,11)==0
replace g = "" if mod(id,31)==0
gen str80 text = cond(mod(id,7)==0,"", cond(mod(id,3)==0,"alpha", ///
    cond(mod(id,3)==1,"  alpha|beta  ","alpha||東京|tail")))
forvalues j=1/8 {
    gen double v`j' = rnormal()
    gen str24 s`j' = "unused payload" + string(mod(id,100))
}
local comparisons=0
foreach selection in "" "if mod(id,3)==0" "in 13/17000" ///
    "if mod(id,3)==1 in 13/17000" "if 0" "if missing(x)" {
    foreach grouping in "" "by(g)" "by(x)" "by(y)" {
        ipolate y x `selection', gen(expected) epolate `grouping'
        cipolate y x `selection', gen(actual) epolate `grouping'
        assert reldif(actual,expected)<1e-12 | (missing(actual)&missing(expected))
        assert id==_n
        drop expected actual
        local ++comparisons
    }
    capture split text `selection', gen(expected) parse("|")
    local expected_rc=_rc
    capture csplit text `selection', gen(actual) parse("|")
    assert _rc==`expected_rc'
    if !`expected_rc' {
        unab expected_vars : expected*
        unab actual_vars : actual*
        local fields : word count `expected_vars'
        assert `fields'==`: word count `actual_vars''
        forvalues j=1/`fields' {
            assert actual`j'==expected`j'
            assert "`: type actual`j''"=="`: type expected`j''"
        }
        drop actual* expected*
    }
    assert id==_n
    local ++comparisons
}
* Empty strings are excluded; whitespace-only strings are valid split inputs.
replace text=""
capture csplit text, gen(actual)
assert _rc==2000
capture confirm variable actual1
assert _rc==111
replace text="   "
split text, gen(expected)
csplit text, gen(actual)
assert actual1==expected1
drop expected* actual*
local comparisons=`comparisons'+2
* Every missing x value must be written even with an uninitialized destination.
replace x=.z
cipolate y x, gen(actual) epolate
assert actual==.
drop actual
local ++comparisons

* Grouped winsor output precision, metadata and missing values match winsor2.
* No-group comparisons use the explicit order-statistic formula, since native
* winsor2 rounds its percentile through a macro (a pre-existing difference).
gen byte vb=mod(id,103)
gen int vi=mod(id*37,20000)
replace vb=.z in 7
replace vi=.a in 9
foreach storage in float double {
    set type `storage'
    foreach selection in "" "if mod(id,3)==0" "in 13/17000" "if id>50 in 13/17000" {
        replace v1=.a in 3
        replace v2=.z in 4
        label variable v1 "first source"
        format v1 %15.4f
        cwinsor v1 v2 vb vi, by(g) suffix(_all)
        winsor2 v1 v2 vb vi, by(g) suffix(_ref)
        foreach v in v1 v2 vb vi {
            assert `v'_all==`v'_ref
            drop `v'_all `v'_ref
        }
        cwinsor v1 v2 vb vi `selection', by(g) suffix(_actual)
        winsor2 v1 v2 vb vi `selection', by(g) suffix(_expected)
        foreach v in v1 v2 vb vi {
            assert `v'_actual==`v'_expected
            assert "`: type `v'_actual'"=="`: type `v'_expected'"
            assert "`: format `v'_actual'"=="`: format `v'_expected'"
            assert `"`: variable label `v'_actual'"'==`"`: variable label `v'_expected'"'
        }
        drop v1_actual v2_actual vb_actual vi_actual v1_expected v2_expected vb_expected vi_expected
        local comparisons=`comparisons'+2
    }
}
set type float
* All outputs belong to the wrapper transaction: direct-store errors must undo them.
foreach fault in outputs store post {
    quietly datasignature
    local before `r(datasignature)'
    global CTOOLS_cwinsor_fault `fault'
    capture cwinsor v1 v2, by(g) suffix(_failure)
    local rc=_rc
    global CTOOLS_cwinsor_fault
    assert `rc'==cond("`fault'"=="store",920,1)
    quietly datasignature
    assert "`before'"=="`r(datasignature)'"
    local ++comparisons
}
gen float vf=v2
foreach selection in "" "if mod(id,3)==0" "in 13/17000" {
    gen byte selected=0
    replace selected=1 `selection'
    foreach source in v1 vf {
        _pctile `source' if selected, p(5 95)
        * winsor2 expands numeric locals, rounding the ungrouped cutoffs.
        local lower=r(r1)
        local upper=r(r2)
        local storage : type `source'
        gen `storage' expected=min(max(`source',`lower'),`upper') ///
            if selected & !missing(`source')
        cwinsor `source' `selection', cuts(5 95) suffix(_direct)
        assert `source'_direct==expected
        drop `source'_direct expected
        local ++comparisons
    }
    drop selected
}
* Missing-x rows must be excluded before attempting SPI reads of long by-keys.
clear
set obs 3
gen double x=_n
gen double y=2*x
replace x=.z in 2
gen strL longkey="g"
forvalues j=1/12 {
    quietly replace longkey=longkey+longkey in 2
}
cipolate y x, by(longkey) gen(actual) epolate
assert actual==y if !missing(x)
assert actual==. if missing(x)
drop actual
replace x=.
cipolate y x, by(longkey) gen(actual)
assert actual==.
local comparisons=`comparisons'+2
di "PERFORMANCE_PATHS_COMPARISONS=`comparisons'"
di "PERFORMANCE_PATHS_COMPLETE"
    global CTOOLS_CASE_COMPLETE "performance_paths"
    exit
}

if "`__case'" == "refactor_timing" {
adopath ++ build
clear
set obs 200
gen double x = mod(_n, 17)
gen double y = _n
gen int g = mod(_n, 4)
timer clear 99
timer on 99
quietly summarize y
local mean = r(mean)
quietly _ctools_clock now
assert r(mean) == `mean'
assert `now' > 0
csort x, verbose
creghdfe y x, absorb(g) verbose
timer off 99
quietly timer list 99
assert r(nt99) == 1
assert r(t99) > 0
local prior = r(t99)
csort y, verbose
quietly timer list 99
assert r(nt99) == 1
assert r(t99) == `prior'
di "CTOOLS_REFACTOR_TIMING_COMPLETE"
    global CTOOLS_CASE_COMPLETE "refactor_timing"
    exit
}

if "`__case'" == "threads" {
/* Thread settings: ctools, threads() session and saved settings, and how the
   dispatcher applies them. The developer's own session and saved settings are
   restored at the end. */
do validation/support/validate_setup.do
clear
set obs 2000
set seed 4412
gen double x = runiform()

* csort skips the plugin when x is already sorted, so reshuffle before each call.
capture program drop shuffle_x
program define shuffle_x
    quietly replace x = runiform()
end

local orig_session "$CTOOLS_THREADS"
global CTOOLS_THREADS
quietly ctools, environment_check
local orig_saved "`r(threads_saved)'"
local limit = r(threads_limit)
local cap = r(threads_cap)
local env_threads : env CTOOLS_NUM_THREADS

* Defaults: capped, never above the usable CPUs.
if `cap' == 16 & r(threads) >= 1 & r(threads) <= `limit' & r(cpus_available) <= r(cpus_detected) {
    test_pass "default is capped and within usable CPUs"
}
else test_fail "default thread policy" "threads=`r(threads)' limit=`limit' cap=`cap'"
if "`r(threads_source)'" != "session setting" & r(threads) == r(threads_default) {
    test_pass "no session setting after clearing"
}
else test_fail "cleared session setting" "source=`r(threads_source)'"

* Session setting applies to commands without threads().
quietly ctools, threads(2)
local expect = min(2, `limit')
if "$CTOOLS_THREADS" == "2" & r(threads) == `expect' & "`r(threads_source)'" == "session setting" {
    test_pass "session setting stored"
}
else test_fail "session setting" "global=$CTOOLS_THREADS threads=`r(threads)'"
shuffle_x
csort x
if r(threads_max) == `expect' test_pass "commands use the session setting"
else test_fail "session setting reaches commands" "threads_max=`r(threads_max)'"
shuffle_x
csort x, threads(1)
if r(threads_max) == 1 test_pass "threads() overrides the session setting"
else test_fail "threads() override" "threads_max=`r(threads_max)'"
shuffle_x
csort x, threads(256)
if r(threads_max) == `limit' test_pass "explicit request clamped to usable CPUs"
else test_fail "explicit clamp" "threads_max=`r(threads_max)' limit=`limit'"

* Invalid settings fail clearly and can be cleared.
global CTOOLS_THREADS abc
shuffle_x
capture csort x
if _rc == 198 test_pass "invalid session global rejected"
else test_fail "invalid session global" "rc=`=_rc'"
quietly ctools, environment_check
if "`r(threads_session)'" == "invalid" test_pass "invalid session global reported"
else test_fail "invalid session report" "session=`r(threads_session)'"
quietly ctools, threads(default)
if "$CTOOLS_THREADS" == "" & "`r(threads_source)'" != "session setting" test_pass "threads(default) clears the session"
else test_fail "threads(default)" "global=$CTOOLS_THREADS"
shuffle_x
csort x
test_pass "commands run after clearing"
foreach bad in 0 257 -1 1.5 abc "2 3" {
    capture ctools, threads(`bad')
    if _rc == 198 test_pass "rejected threads(`bad')"
    else test_fail "threads(`bad')" "rc=`=_rc'"
}
capture ctools, permanently
if _rc == 198 test_pass "permanently requires threads()"
else test_fail "permanently alone" "rc=`=_rc'"

* Saved setting: written outside the ado tree, applies in this session once
* the session setting is cleared, and is removed by threads(default).
if "`orig_saved'" == "invalid" {
    test_skip "saved setting" "existing settings file has an invalid threads= value"
}
else {
    quietly ctools, threads(3) permanently
    local file `"`r(threads_file)'"'
    capture confirm file `"`file'"'
    if !_rc & "`r(threads_saved)'" == "3" & "$CTOOLS_THREADS" == "3" test_pass "saved setting written"
    else test_fail "saved setting written" "saved=`r(threads_saved)' file=`file'"
    quietly ctools, threads(default)
    if "`env_threads'" == "" {
        if r(threads) == min(3, `limit') & "`r(threads_source)'" == "saved setting" {
            test_pass "saved setting is the default"
        }
        else test_fail "saved default" "threads=`r(threads)' source=`r(threads_source)'"
    }
    quietly ctools, threads(default) permanently
    if "`r(threads_saved)'" == "" & "$CTOOLS_THREADS" == "" test_pass "saved setting removed"
    else test_fail "saved setting removed" "saved=`r(threads_saved)'"
    if "`orig_saved'" != "" quietly ctools, threads(`orig_saved') permanently
}

global CTOOLS_THREADS `orig_session'
print_summary "threads"
global CTOOLS_COMPONENT_COMPLETE "threads"
    global CTOOLS_CASE_COMPLETE "threads"
    exit
}

if "`__case'" == "transport_strl" {
/* Real-SPI strL regression. Supply the standalone benchmark_transport plugin.
   The driver must invoke this through the stata shell alias and finish with exit, clear. */
args transport_plugin
clear all
set more off
program transport_io, plugin using("`transport_plugin'")
forvalues threads = 1(7)8 {
    foreach width in 8 244 2045 {
        clear
        quietly set obs 25000
        quietly gen strL text1 = "x"
        forvalues power = 1/11 {
            quietly replace text1 = text1 + text1
        }
        quietly replace text1 = substr(string(_n,"%06.0f") + text1,1,`width')
        quietly gen double number = cond(mod(_n,19)==0,.z,_n/7)
        quietly gen strL text2 = "café_" + string(_n)
        quietly replace text1 = "" if mod(_n,17)==0
        quietly gen str12 fixed = string(_n)
        quietly datasignature
        local signature "`r(datasignature)'"
        foreach vars in "text1" "text1 text2" "number text1 fixed text2" "text1 number text2 fixed" {
            foreach rep in 1 2 3 {
                plugin call transport_io `vars', "strl" "`threads'" "read" "1"
                plugin call transport_io `vars' if mod(_n,3)!=0 in 3/24999, "strl" "`threads'" "filteredread" "1"
            }
        }
        quietly datasignature
        assert "`r(datasignature)'" == "`signature'"
    }
}
/* Cross production scheduling thresholds with one to three columns. */
clear
quietly set obs 250003
quietly gen strL text = "café_" + string(_n)
quietly gen strL text2 = "second_" + string(_n)
quietly gen double number = cond(mod(_n,19)==0,.z,_n/7)
foreach threads in 1 8 12 {
    foreach vars in "text" "text text2" "text text2 number" {
        local __ctools_strw "64,64,0"
        plugin call transport_io `vars', "large_strl" "`threads'" "read" "1"
        plugin call transport_io `vars' if mod(_n,3)!=0, "large_strl" "`threads'" "filteredread" "1"
        local __ctools_strw ""
        plugin call transport_io `vars', "large_strl" "`threads'" "read" "1"
    }
}
/* Bounds and write rejection must remain explicit errors. */
clear
quietly set obs 3
quietly gen strL text = "x"
forvalues power = 1/11 {
    quietly replace text = text + text
}
quietly replace text = substr(text,1,2046)
capture plugin call transport_io text, "oversize" "8" "read" "0"
assert _rc == 5
quietly replace text = "short"
capture plugin call transport_io text, "write" "8" "identity" "0"
assert _rc == 5
display "TRANSPORT_STRL_COMPLETE"
    global CTOOLS_CASE_COMPLETE "transport_strl"
    exit
}

di as error "Unknown contracts case: `__case'"
exit 198
