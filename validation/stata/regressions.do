* regressions: named cases share this file, each do-call has its own local scope.
gettoken __case 0 : 0
macro shift
if "`__case'" == "" local __case "regressions"

if "`__case'" == "audit_p1" {
/* Targeted regression tests for audit findings A01--A17. Run from the root. */
do validation/support/validate_setup.do

program p1_helpers
    capture noisily do validation/stata/oracle_tests.do validation_helpers
    local rc = _rc
    * Even a failed mutation assertion must not leak stubs into later cases.
    foreach cmd in csort cimport qreg cqreg ppmlhdfe cppmlhdfe crangestat rangestat helper_ppml_result helper_expect {
        capture program drop `cmd'
    }
    macro drop HELPER_MUTATION
    exit `rc'
end

program p1_merge
    clear
    set obs 2
    gen id = _n
    gen val = 100 + _n
    tempfile usingdata
    save `usingdata'
    replace val = 999
    tempname master sentinel
    frame copy default `master'
    frame create `sentinel'
    frame change `master'
    drop val
    cmerge 1:1 id using `usingdata'
    assert "`c(frame)'" == "`master'"
    assert val == 100 + id
    frame default: assert val == 999
    drop _merge
    capture cmerge 1:1 id using `usingdata', keepusing(not_a_variable)
    assert _rc != 0
    assert "`c(frame)'" == "`master'"
    frame default: assert val == 999
    frame change default
end

program p1_sample
    foreach qualifier in "if id<=5" "in 1/5" "if id<=7 in 1/5" {
        foreach pct in 0 100 {
            clear
            set obs 10
            gen id = _n
            gen __csample_keep__ = 123
            csample `pct' `qualifier'
            assert __csample_keep__ == 123
            count if id > 5
            assert r(N) == 5
            assert _N == 5 + (`pct' == 100)*5
        }
    }
end

program p1_groups
    foreach cmd in csample cbsample {
        foreach kind in string numeric mixed missing {
            clear
            set obs 10
            gen id = _n
            gen str1 s = cond(_n<=5,"A","B")
            gen double g = cond(_n<=5,.1,1e12)
            if "`kind'" == "missing" replace g = cond(_n<=5,.a,.b)
            local groups s
            if inlist("`kind'","numeric","missing") local groups g
            if "`kind'" == "mixed" local groups s g
            if "`cmd'" == "cbsample" & "`kind'" == "missing" {
                replace g = 1 in 1/5
                cbsample 1, strata(g)
                assert _N == 1 & g == 1
                continue
            }
            if "`cmd'" == "csample" csample, count(1) by(`groups')
            else cbsample 1, strata(`groups')
            assert _N == 2
            isid `groups', missok
        }
    }
end

program p1_matching
    foreach option in "" "noreplacement" "noreplacement descending" {
        clear
        set obs 3
        gen byte treat = _n==3
        gen double ps = cond(_n==1,.1,cond(_n==2,.51,.5))
        gen double y = cond(_n==1,10,cond(_n==2,20,100))
        cpsmatch treat, pscore(ps) outcome(y) `option'
        assert r(att) == 80
        assert _weight[2] == 1 & _weight[3] == 1
        assert _support == 1
        assert _id[3] == 2
    }
    clear
    set obs 4
    gen byte treat = _n>=3
    gen double ps = cond(_n==1,.1,cond(_n==2,.51,cond(_n==3,.11,.5)))
    gen double y = cond(_n==1,10,cond(_n==2,20,cond(_n==3,100,200)))
    cpsmatch treat, pscore(ps) outcome(y) noreplacement
    assert _id[3] == 1 & _id[4] == 2
    assert _weight == 1
    assert r(att) == 135
    cpsmatch treat, pscore(ps) outcome(y) noreplacement caliper(.001)
    assert _support == 0 if treat
    capture cpsmatch treat, pscore(ps) noreplacement neighbor(2)
    assert _rc == 198
    capture cpsmatch treat, pscore(ps) noreplacement kernel
    assert _rc == 198
    capture cpsmatch treat, pscore(ps) noreplacement radius
    assert _rc == 198
    clear
    set obs 5
    gen byte treat = _n>2
    gen double ps = cond(_n==1,.25,cond(_n==2,.75,.5))
    cpsmatch treat, pscore(ps) noreplacement
    assert _id[3]==1 & _id[4]==2
    assert _support[5]==0
    assert r(n_matched)==2
    clear
    set seed 12345
    set obs 500
    gen x1 = rnormal()
    gen x2 = rnormal()
    gen byte treat = runiform() < normal(.3*x1+.2*x2)
    gen double y = 2*treat+x1+x2+rnormal()
    gen double saved_y = y
    cpsmatch treat x1 x2, outcome(y)
    assert y == saved_y
    assert inrange(_pscore,0,1)
    assert _support == 1
    assert !missing(_weight, r(att))
end

program p1_labels
    clear
    set obs 3
    gen str1 s = char(64+_n)
    label define existing 10 "A" 20 "B" 30 "unused"
    gen other = 30
    label values other existing
    encode s, gen(expected) label(existing)
    cencode s, gen(actual) label(existing)
    assert actual == expected
    local unused : label existing 30
    assert "`unused'" == "unused"
    cencode s, gen(again) label(existing)
    assert again == actual
    gen str1 s2 = s
    cencode s s2, gen(multi1 multi2) label(existing)
    assert multi1 == expected & multi2 == expected
end

program p1_bins
    clear
    set obs 40
    gen double x = _n
    gen double y = 2*x + mod(_n,3)
    gen double g = cond(_n<=20,10,20)
    cbinscatter y x, by(g) nquantiles(4) nograph
    assert e(N) == 40 & e(num_groups) == 2
    matrix bins = e(bindata)
    * binscatter bins every group on the pooled quartiles of x (10, 20, 30, 40):
    * the groups' disjoint x ranges each fill two of them, and the remaining
    * rows of the 4-bins-per-group matrix stay missing
    assert rowsof(bins) == 8
    scalar total = 0
    scalar filled = 0
    forvalues i=1/8 {
        if missing(bins[`i',1]) continue
        scalar filled = filled + 1
        assert !missing(bins[`i',3],bins[`i',4],bins[`i',5])
        assert bins[`i',1] == 1 + (`i' > 2) & bins[`i',5] == 10
        scalar total = total + bins[`i',5]
    }
    assert filled == 4
    assert total == 40
    replace g = cond(_n<=20,-.1,1e12)
    cbinscatter y x, by(g) nquantiles(4) nodraw
    assert e(num_groups) == 2
    gen str1 gs = cond(_n<=20,"A","B")
    cbinscatter y x, by(gs) nquantiles(4) nograph
    assert e(num_groups) == 2
end

program p1_import
    tempfile csv xlsx
    foreach encoding in ascii utf8 {
        foreach width in 2044 2045 2046 3000 {
            foreach format in csv xlsx {
                clear
                set obs 1
                gen strL text = ""
                if "`encoding'" == "ascii" {
                    mata: st_sstore(1,"text",invtokens(J(1,`width',"a"),""))
                }
                else {
                    local symbol = uchar(233)
                    local chars = floor(`width'/2)
                    mata: st_sstore(1,"text",invtokens(J(1,`chars',"`symbol'"),""))
                    if mod(`width',2) replace text = text + "a"
                }
                local expected = text[1]
                if "`format'" == "csv" export delimited using `csv', replace
                else export excel using "`xlsx'.xlsx", firstrow(variables) replace
                clear
                set obs 1
                gen sentinel = 42
                if "`format'" == "csv" {
                    capture noisily cimport delimited using `csv', varnames(1) clear
                    local rc = _rc
                }
                else {
                    capture noisily cimport excel using "`xlsx'.xlsx", firstrow clear
                    local rc = _rc
                }
                * CSV header inference is ambiguous for a single string column;
                * explicitly request its header. Both readers now support strL.
                assert `rc' == 0
                assert _N == 1
                assert strlen(text) == `width'
                assert text == "`expected'"
                if `width' > 2045 {
                    local storage : type text
                    assert "`storage'" == "strL"
                }
            }
        }
    }
    erase "`xlsx'.xlsx"
end

program p1_iv
    clear
    set seed 2050
    set obs 400
    gen group = ceil(_n/20)
    replace group = 999 in 400
    gen double z = rnormal()
    gen double x = z + rnormal()
    gen double y = 10*group + 2*x + rnormal()
    replace z = . in 6
    civreghdfe y (x=z) if _n>3, absorb(group) residuals(savedres)
    count if e(sample)
    assert r(N) == e(N) & r(N) == 395
    assert !missing(savedres) if _n>3 & _n!=6 & _n!=400
    assert missing(savedres) if _n<=3 | _n==6 | _n==400
    assert !e(sample) if _n<=3 | _n==6 | _n==400
    assert "`e(resid)'" == "savedres"
    predict double residual, residuals
    predict double fit, xbd
    assert residual == savedres
    assert abs(y-fit-residual) < 1e-10 if e(sample)
    tempvar meanres
    egen double `meanres' = mean(residual), by(group)
    assert abs(`meanres') < 1e-7 if e(sample)
    estimates store iv_saved
    quietly regress y x
    estimates restore iv_saved
    predict double restored, residuals
    assert restored == savedres
    gen _civreghdfe_resid = 999
    civreghdfe y (x=z), absorb(group) residuals2
    count if e(sample)
    assert r(N) == e(N) & r(N) == 398
    assert !missing(_civreghdfe_resid) if e(sample)
    assert _civreghdfe_resid != 999 if e(sample)
end

program p1_variance
    foreach n in 10 200 {
        foreach exclude in "" "excludeself" {
            clear
            set obs `n'
            gen time = _n
            gen double x = 1e12 + _n
            gen double small = _n
            crangestat (variance) vx=x (variance) vs=small (sd) sx=x, interval(time . .) `exclude'
            assert abs(vx-vs) < 1e-8
            assert abs(sx^2-vx) < 1e-8
            summarize small
            if "`exclude'" == "" assert abs(vx-r(Var)) < 1e-8
            crangestat (variance) wx=x (variance) ws=small, interval(time -2 2) `exclude'
            assert abs(wx-ws) < 1e-8
        }
    }
    clear
    set obs 400000
    gen time = _n
    gen double x = 1e12 + _n
    crangestat (variance) vx=x, interval(time -99 0)
    assert abs(vx-100*101/12) < 1e-8 in 399900/400000

end

program p1_quantile
    clear
    set seed 54321
    set obs 400
    gen group = ceil(_n/20)
    gen double x = rnormal()
    gen double y = 2*x + rnormal()
    capture cqreg y x, absorb(group)
    assert _rc == 198
    capture cqreg y x, maxiter(1)
    assert _rc == 430
    cqreg y x
    assert e(convcode) == 0
    assert e(iterations) >= 0 & e(iterations) <= 200
    qreg y x
    scalar expected = _b[x]
    cqreg y x
    assert reldif(_b[x],expected) < 1e-7
end

/* The helper self-tests intentionally record failures, then reset the counters. */
clear
set obs 1
gen a = .
gen b = 42
assert_var_equal a b 7 "missing mismatch"
assert $TESTS_FAILED == 1
matrix a = (1,2)
matrix b = (1,2,3)
assert_matrix_equal a b 7 "dimension mismatch"
assert $TESTS_FAILED == 2
matrix_min_sigfigs a b
assert r(min_sigfigs) == 0
capture print_summary "expected helper failures"
assert _rc == 9
global TESTS_PASSED = 0
global TESTS_FAILED = 0
global TESTS_TOTAL = 0
global FAILURE_COUNT = 0

foreach case in helpers merge sample groups matching labels bins import iv variance quantile {
    capture noisily p1_`case'
    if _rc test_fail "P1 `case'" "r(`=_rc')"
    else test_pass "P1 `case'"
}
print_summary "audit P1 regressions"

* Reached only after the complete component script.
global CTOOLS_COMPONENT_COMPLETE "audit_p1"
    global CTOOLS_CASE_COMPLETE "audit_p1"
    exit
}

if "`__case'" == "audit_p2" {
/* Targeted regression tests for audit findings A18--A29. */
do validation/support/validate_setup.do

program p2_labels
    clear
    set obs 8
    gen str2045 source = ""
    global AUDIT_WORD changed
    replace source = char(36)+"AUDIT_WORD" in 1
    replace source = char(96)+"local_name"+char(39) in 2
    replace source = char(96)+char(34)+"quoted"+char(34)+char(39) in 3
    replace source = "a"+char(9)+"b"+char(10)+"c"+char(13)+"d" in 4
    replace source = "café 日本語" in 5
    replace source = char(34)+"plain quotes"+char(34) in 6
    mata: st_sstore(7,"source",invtokens(J(1,2045,"x"),""))
    cencode source, gen(code)
    decode code, gen(roundtrip)
    assert source == roundtrip
    macro drop AUDIT_WORD
    clear
    set obs 4
    gen str1 s = char(64+_n)
    label define signed -10 "A" 0 "B" 7 "C"
    capture encode s, gen(expected) label(signed) noextend
    local native_rc = _rc
    capture cencode s, gen(actual) label(signed) noextend
    assert _rc == `native_rc' & `native_rc' == 459
    capture drop expected actual
    drop in 4
    encode s, gen(expected) label(signed) noextend
    cencode s, gen(actual) label(signed) noextend
    assert actual == expected
end

program p2_clusters
    clear
    set seed 45021
    set obs 400
    gen double x = rnormal()
    gen double y = .7*x+rnormal()
    gen long g = ceil(_n/20)
    gen double fraction = g/100
    gen double large = 1e12+g
    gen str20 textgroup = string(g)
    cqreg y x, vce(cluster g)
    tempname original_b original_V
    matrix `original_b' = e(b)
    matrix `original_V' = e(V)
    foreach group in fraction large textgroup {
        cqreg y x, vce(cluster `group')
        assert e(N_clust) == 20
        assert_matrix_equal `original_b' e(b) $DEFAULT_SIGFIGS "cluster relabeling coefficients"
        assert_matrix_equal `original_V' e(V) $DEFAULT_SIGFIGS "cluster relabeling covariance"
    }
end

program p2_center
    clear
    set obs 30
    gen t = _n
    gen y = sin(t)
    gen x = cos(t)
    gen z = t^2
    tsset t
    capture civreghdfe y (x=z), robust bw(3) center
    assert _rc == 198
end

program p2_import
    clear
    set obs 2
    gen str5 name = cond(_n==1,"Alice","Bob")
    gen value = _n
    tempfile base
    local csv "`base' with spaces.csv"
    local xlsx "`base' with spaces.xlsx"
    export delimited using "`csv'", replace
    export excel using "`xlsx'", firstrow(variables) replace
    foreach syntax in positional using compound {
        if "`syntax'" == "positional" cimport delimited "`csv'", clear
        if "`syntax'" == "using" cimport delimited using "`csv'", clear
        if "`syntax'" == "compound" cimport delimited `"`csv'"', clear
        assert _N == 2 & name[1] == "Alice" & value[2] == 2
        if "`syntax'" == "positional" cimport excel "`xlsx'", firstrow clear
        if "`syntax'" == "using" cimport excel using "`xlsx'", firstrow clear
        if "`syntax'" == "compound" cimport excel `"`xlsx'"', firstrow clear
        assert _N == 2 & name[1] == "Alice" & value[2] == 2
    }
    foreach kind in delimited excel {
        capture cimport `kind' "`csv'" extra, clear
        assert _rc == 198
        assert _N == 2 & name[1] == "Alice"
    }
    foreach endian in le be {
        tempfile utf32
        tempname fh
        file open `fh' using "`utf32'", write binary replace
        if "`endian'" == "le" file write `fh' %1bu (255) %1bu (254) %1bu (0) %1bu (0)
        else file write `fh' %1bu (0) %1bu (0) %1bu (254) %1bu (255)
        foreach byte in 97 44 98 10 49 44 50 10 {
            if "`endian'" == "le" file write `fh' %1bu (`byte') %1bu (0) %1bu (0) %1bu (0)
            else file write `fh' %1bu (0) %1bu (0) %1bu (0) %1bu (`byte')
        }
        file close `fh'
        * UTF-32 is supported: verify decoded contents for automatic and
        * explicit endian selection, rather than expecting obsolete errors.
        foreach enc in auto utf-32 utf-32`endian' {
            local encopt
            if "`enc'" != "auto" local encopt encoding(`enc')
            cimport delimited using "`utf32'", clear `encopt'
            assert _N == 1 & c(k) == 2 & a == 1 & b == 2
        }
    }
    capture cimport delimited "`csv'", clear encoding(not-a-charset)
    assert _rc == 198
    assert _N == 1 & c(k) == 2 & a == 1 & b == 2
    erase "`csv'"
    erase "`xlsx'"
end

program p2_sort
    clear
    set obs 600
    gen id = _n
    gen key = mod(600-_n,17)
    gen double nonkey = _n^2
    gen str20 text = "row"+string(_n)
    foreach qualifier in "if id<=2" "in 1/2" {
        capture csort key `qualifier'
        assert _rc == 101
        assert id == _n & nonkey == _n^2
    }
    tempfile source expected
    save `source'
    sort key, stable
    gen rank = _n
    save `expected'
    foreach algorithm in auto lsd msd merge timsort sample counting ips4o {
        use `source', clear
        csort key, algorithm(`algorithm') nosortedby
        gen rank = _n
        cf _all using `expected'
    }
    use `source', clear
    csort key, stream(1) nosortedby
    gen rank = _n
    cf _all using `expected'
end

program p2_sample_by
    foreach pct in 30 100 {
        clear
        set obs 100
        gen long id = _n
        gen byte g = ceil(_n/10)
        tempfile data expected
        save `data'
        sample `pct', by(g)
        collapse (count) n=id, by(g)
        save `expected'
        use `data', clear
        csample `pct', by(g)
        collapse (count) n=id, by(g)
        cf _all using `expected'
    }
end

program p2_io
    clear
    set obs 3
    gen strL text = "long string"
    gen id = _n
    cencode text, gen(encoded)
    decode encoded, gen(decoded)
    assert decoded == text
    drop encoded decoded
    label drop encoded
    capture csort id, stream(1) nosortedby
    assert _rc != 0
    assert text == "long string" & id == _n
    capture csort id, nosortedby
    assert _rc != 0
    assert text == "long string" & id == _n
    mata: st_sstore(., "text", J(3,1,invtokens(J(1,3000,"x"),"")))
    capture cencode text, gen(encoded)
    assert _rc != 0
    assert strlen(text) == 3000 & id == _n
end

program p2_identity
    ctools, version
    assert "`r(version)'" != ""
    assert "`r(version)'" == "`r(plugin_version)'"
    assert "`r(plugin_revision)'" != ""
end

foreach case in labels clusters center import sort sample_by io identity {
    capture noisily p2_`case'
    if _rc test_fail "A18-A29: `case'" "r(`=_rc')"
    else test_pass "A18-A29: `case'"
}
print_summary "P2 audit regressions"

* Reached only after the complete component script.
global CTOOLS_COMPONENT_COMPLETE "audit_p2"
    global CTOOLS_CASE_COMPLETE "audit_p2"
    exit
}

if "`__case'" == "oct21" {
version 14.1
args fixtures plugin_dir
adopath ++ build
if `"`plugin_dir'"' != "" adopath ++ `"`plugin_dir'"'

* CRC and XML corruption fail before replacing the caller's data.
foreach name in bitflip truncated-xml {
    clear
    set obs 1
    gen sentinel = 42
    capture noisily cimport excel `"`fixtures'/`name'.xlsx"', clear
    assert _rc != 0
    assert _N == 1 & c(k) == 1
    assert sentinel == 42
}
cimport excel `"`fixtures'/valid.xlsx"', clear
assert _N == 2 & A[1] == 123 & A[2] == 456

* An existing companion must not leave a newly published geometry file.
cimport shp `"`fixtures'/point.shp"', clear
capture noisily cexport shp using `"`fixtures'/partial.shp"', shx
assert _rc == 602
capture confirm file `"`fixtures'/partial.shp"'
assert _rc == 601
cexport shp using `"`fixtures'/complete.shp"', shx
confirm file `"`fixtures'/complete.shx"'

* Explicit protection precedes scans for every binary reader.
foreach kind in dbase shp sas spss sasxport5 sasxport8 {
    clear
    set obs 2
    gen sentinel = 43 - _n
    sort sentinel
    label var sentinel "snapshot label"
    char _dta[protected] "yes"
    capture noisily cimport `kind' `"`fixtures'/corrupt.dbf"', clear snapshot
    assert _rc != 0
    assert _N == 2 & c(k) == 1
    assert sentinel[1] == 41 & sentinel[2] == 42
    local label : variable label sentinel
    assert `"`label'"' == "snapshot label"
    local protected : char _dta[protected]
    assert "`protected'" == "yes"
}
* Preserve/restore on failures after the scan as well as before it.
clear
set obs 2
gen x = _n
cexport dbase using `"`fixtures'/valid.dbf"'
replace x = 42
global CTOOLS_cimport_fault load
capture noisily cimport dbase `"`fixtures'/valid.dbf"', clear snapshot
local failed = _rc
macro drop CTOOLS_cimport_fault
assert `failed' != 0
assert _N == 2 & x[1] == 42 & x[2] == 42
cimport dbase `"`fixtures'/valid.dbf"', clear snapshot
assert x[1] == 1 & x[2] == 2

* A blocked second companion must leave neither the transport nor the first companion.
format x %9.0g
label define yesno 1 "one" 2 "two"
label values x yesno
capture noisily cexport sasxport5 using `"`fixtures'/transport.xpt"', vallabfile(both)
assert _rc == 602
foreach file in transport.xpt formats.xpf {
    capture confirm file `"`fixtures'/`file'"'
    assert _rc == 601
}
label values x

* Unicode paths exercise staging, direct IO, vendor codecs and commit.
local unicode = uchar(233) + uchar(28204) + uchar(35430)
mkdir `"`fixtures'/`unicode'"'
foreach format in delimited dbase excel sas spss sasxport8 {
    local ext csv
    if "`format'" == "dbase" local ext dbf
    if "`format'" == "excel" local ext xlsx
    if "`format'" == "sas" local ext sas7bdat
    if "`format'" == "spss" local ext sav
    if "`format'" == "sasxport8" local ext v8xpt
    local filename `"`fixtures'/`unicode'.`ext'"'
    * XPORT derives an ASCII SAS member name from the basename.
    if "`format'" == "sasxport8" local filename `"`fixtures'/`unicode'/data.`ext'"'
    cexport `format' using `"`filename'"', replace
    cimport `format' `"`filename'"', clear
    assert _N == 2
}

cexport excel using `"`fixtures'/`unicode'.xls"', replace
cimport excel `"`fixtures'/`unicode'.xls"', clear
assert _N == 2

* Unicode whitespace matches native destring without discarding interior text.
clear
set obs 10
gen str30 s = ""
local i 0
foreach cp in 9 160 5760 8192 8195 8201 8232 8239 8287 12288 {
    local ++i
    replace s = uchar(`cp') + "123" + uchar(`cp') in `i'
}
destring s, gen(native) force
cdestring s, gen(ct) force
assert native == ct

* Every affected prediction path honors if/in and publishes only on success.
clear
set obs 100
gen double x = _n/100
gen double y = mod(_n,4)
cppmlhdfe y x, d(dsum)
foreach opt in mu eta xbd response scores pearson deviance anscombe working xb d stdp {
    predict double all_`opt', `opt'
    predict double sub_`opt' if _n < 3 in 1/5, `opt'
    assert all_`opt' == sub_`opt' in 1/2
    assert missing(sub_`opt') in 3/100
}
drop dsum
capture noisily predict double unfinished if _n < 3, mu
assert _rc != 0
capture confirm variable unfinished
assert _rc != 0

di "CTOOLS_OCT21_COMPLETE"
    global CTOOLS_CASE_COMPLETE "oct21"
    exit
}

if "`__case'" == "oct22" {
* Regression coverage for audit A01-A29. Run from the repository root.
version 16
clear all
set more off
set seed 9854
adopath ++ build
local candidate : environment CTOOLS_TEST_PLUGIN_DIR
if `"`candidate'"' != "" adopath ++ `"`candidate'"'
args fixtures
if `"`fixtures'"' == "" local fixtures "validation/fixtures/oct22"

* A01, A11-A14: independent IV reference and positive-DF tiny model.
set obs 100
gen double z=rnormal()
gen double x=z+rnormal()
gen double y=10+.05*x+rnormal()
foreach options in "small" "small robust" "noconstant small" {
    quietly ivregress 2sls y (x=z), `options'
    matrix expected=e(V)
    quietly test x
    scalar expected_F=r(F)
    quietly summarize y
    scalar expected_tss=(r(N)-1)*r(Var)
    if strpos("`options'","noconstant") scalar expected_tss=expected_tss+r(N)*r(mean)^2
    quietly civreghdfe y (x=z), `options'
    matrix actual=e(V)
    noi di "IV `options' Vdiff=" mreldif(actual,expected) " F=" e(F) " expected=" expected_F " TSSdiff=" reldif(e(tss),expected_tss)
    assert mreldif(actual,expected)<1e-8
    assert reldif(e(F),expected_F)<1e-8
    assert reldif(e(tss),expected_tss)<1e-12
    assert e(df_m)==1
}
quietly civreghdfe y (x=z)
matrix large=e(V)
scalar rms=e(rmse)
quietly civreghdfe y (x=z), small
matrix small=e(V)
assert reldif(small[1,1]/large[1,1],100/98)<1e-10
quietly civreghdfe y (x=z), small dofminus(2) sdofminus(3)
matrix adjusted=e(V)
assert e(df_r)==93
assert reldif(adjusted[1,1]/small[1,1],98/93)<1e-10
quietly civreghdfe y (x=z), robust
matrix large=e(V)
quietly civreghdfe y (x=z), robust dofminus(2)
matrix adjusted=e(V)
assert reldif(adjusted[1,1]/large[1,1],100/98)<1e-10
clear
input double(y x z)
1 1 1
4 3 2
2 2 3
5 4 4
end
quietly civreghdfe y (x=z), small
assert e(N)==4 & e(df_r)==2

* A02/A16: fitted and kernel cluster covariance, omitted/base rank.
clear
set obs 2000
set seed 9854
gen double x=5*runiform()
gen double y=2+.5*x+(1+4*x)*rnormal()
gen cl=ceil(_n/20)
foreach density in fitted kernel(gaussian) {
    quietly qreg y x, vce(excluster cl, `density')
    matrix expected=e(V)
    quietly cqreg y x, vce(cluster cl) denmethod(`density')
    matrix actual=e(V)
    assert mreldif(actual,expected)<1e-7
}
gen double duplicate=x
gen category=mod(_n,3)
foreach rhs in "x duplicate" "x i.category" {
    quietly qreg y `rhs'
    scalar expected_rank=e(rank)
    quietly cqreg y `rhs'
    assert e(rank)==expected_rank
}
replace cl=1
quietly cqreg y x, vce(cluster cl)
matrix actual=e(V)
scalar undefined_v=actual[1,1]
assert missing(undefined_v)

* A03/A04: invalid metadata/indexes and record-boundary SHP truncation.
foreach bad in truncated_sst truncated_styles truncated_book truncated_rels invalid_index {
    clear
    set obs 1
    gen sentinel=99
    capture cimport excel "`fixtures'/`bad'.xlsx", clear
    assert _rc!=0
    assert sentinel==99 & _N==1
}
capture cimport shp "`fixtures'/truncated.shp", clear
assert _rc!=0
assert _N==0
cimport shp "`fixtures'/good.shp", clear
assert _N==2

* A08/A09: XML text nodes and OOXML strings, native interchange both ways.
foreach textcase in cdata comment shared_cdata shared_comment escaped {
    cimport excel "`fixtures'/`textcase'.xlsx", clear
    if inlist("`textcase'","cdata","shared_cdata") assert A=="customer & name"
    if inlist("`textcase'","comment","shared_comment") assert A=="beforeafter"
    if "`textcase'"=="escaped" assert A=="a_x000D_b"
}
clear
set obs 4
gen str80 text=cond(_n==1,"a"+char(1)+"b", cond(_n==2,"a"+char(13)+"b",cond(_n==3,"_x0041_","_x005F_x0041_")))
tempfile exported native original
save `original'
cexport excel using "`exported'.xlsx", replace
import excel using "`exported'.xlsx", clear
rename A text
cf _all using `original'
use `original', clear
export excel using "`native'.xlsx", replace
cimport excel using "`native'.xlsx", clear
rename A text
cf _all using `original'

* The existing workbook's shared-string indices survive appends/updates.
use `original', clear
cexport excel using "`exported'.xlsx", sheet(Sheet1, modify) cell(B1)
import excel using "`exported'.xlsx", clear
assert A==B
rename A text
drop B
cf _all using `original'
erase "`exported'.xlsx"
erase "`native'.xlsx"

* A05: weight() inside a quoted condition must remain data.
clear
set obs 6
gen str30 s=cond(_n<=3,"weight(w)","other")
gen w=42
cbsample if s=="weight(w)"
assert _N==3 & w==42 & s=="weight(w)"
cbsample, weight(w)
quietly summarize w
assert r(sum)==3

* A06/A07: Unicode codepoint removal and metadata for generate/replace.
clear
set obs 2
gen str20 s=cond(_n==1,"12©","12é")
clonevar expected=s
clonevar actual=s
quietly destring expected, ignore("é¢") replace
quietly cdestring actual, ignore("é¢") replace
assert actual==expected
replace s="12é"
label var s "Original label"
char s[note1] "Original characteristic"
quietly destring s, gen(expected_num) ignore("é¢")
quietly cdestring s, gen(actual_num) ignore("é¢")
assert actual_num==expected_num
assert `"`: variable label actual_num'"'==`"`: variable label expected_num'"'
assert `"`: char actual_num[note1]'"'==`"`: char expected_num[note1]'"'
quietly cdestring s, replace ignore("é¢")
assert s==12
assert `"`: variable label s'"'=="Original label"
assert `"`: char s[note1]'"'=="Original characteristic"

* A10: prefix threshold, multi-row windows and excludeself cancellation.
clear
set obs 64
gen double x=_n-1
replace x=1e16 in 1
gen key=_n
crangestat (mean) m=x (sum) s=x, interval(key 0 0)
assert m==x & s==x
crangestat (mean) m2=x (sum) s2=x, interval(key -1 1) excludeself
assert m2==x & s2==2*x if inrange(_n,3,63)

* A15: every outcome-dependent PPML prediction accepts TS expressions.
clear
set obs 100
gen time=_n
tsset time
gen double x=rnormal()
gen double y=exp(.2*x+rnormal())
quietly cppmlhdfe L.y x, noabsorb
predict double mu, mu
foreach option in response scores pearson deviance anscombe working {
    predict double result, `option'
    assert !missing(result) if e(sample)
    if "`option'"=="response" assert reldif(result,L.y-mu)<1e-12 if e(sample)
    drop result
}

* A17-A19: expanded using varlists, empty _n merge and assert before keep.
clear
set obs 3
gen id=_n
gen a1=10*_n
gen a2=20*_n
gen extra=99
tempfile usingfile
save `usingfile'
foreach keepvars in "a*" "a1-a2" {
    clear
    set obs 3
    gen id=_n
    cmerge 1:1 id using `usingfile', keepusing(`keepvars')
    assert a1==10*id & a2==20*id
    capture confirm variable extra
    assert _rc==111
}
clear
gen master=.
cmerge 1:1 _n using `usingfile', keepusing(a1)
assert _N==3
capture confirm variable extra
assert _rc==111
capture confirm variable id
assert _rc==111
clear
set obs 0
gen id=.
capture cmerge 1:1 id using `usingfile', assert(match) keep(match)
assert _rc==9
assert _N==3 & _merge==2

* A20/A21/A27-A29: user frame ownership, singleton points, graph and metadata.
clear
set obs 8
gen x=mod(_n,2)
gen y=_n
gen str1 g=cond(_n==8,"b","a")
gen category=mod(_n,3)
frame create _ctools_graph_frame
frame _ctools_graph_frame: set obs 1
frame _ctools_graph_frame: gen protected=99
cbinscatter y x, discrete by(g) linetype(none) nograph
frame _ctools_graph_frame: assert protected==99
matrix bins=e(bindata)
matrix list bins
assert rowsof(bins)>=3
assert bins[3,4]==8 & bins[3,5]==1
assert "`e(by)'"=="g"
foreach line in lfit qfit {
    capture cbinscatter y x, discrete by(g) linetype(`line') nograph
    assert _rc==2001
}
cbinscatter y x, controls(i.category) nograph
assert "`e(controls)'"=="i.category"
program define twoway
    assert strpos(`"`0'"',"connect(l)")>0
end
cbinscatter y x, linetype(connect)
cbinscatter y x, by(g) linetype(connect)
program drop twoway
foreach line in none connect {
    capture cbinscatter y x, linetype(`line') reportreg nograph
    assert _rc==198
}
cbinscatter y x, reportreg nograph
assert "`e(cmd)'"=="cbinscatter"
frame drop _ctools_graph_frame

* Reporting must also work with residualization, groups, weights and polynomials.
clear
set obs 100
gen double x=rnormal()
gen double c=rnormal()
gen double y=2*x+c+rnormal()
gen group=mod(_n,2)
gen fe=mod(_n,10)
gen w=1+mod(_n,3)
foreach options in "controls(c)" "controls(c) absorb(fe)" ///
    "controls(c) by(group) linetype(qfit)" "method(binsreg) controls(c) absorb(fe)" {
    quietly cbinscatter y x [aw=w], `options' nograph
    matrix expected_bins=e(bindata)
    cbinscatter y x [aw=w], `options' reportreg nograph
    matrix actual_bins=e(bindata)
    assert mreldif(actual_bins,expected_bins)==0
    assert "`e(cmd)'"=="cbinscatter"
}

* A23-A25: preserve collided inputs; reference domains and treatment coding.
clear
set obs 10
gen t=mod(_n,2)
gen double p=_n/11
gen y=_n
gen double _pscore=p
capture cpsmatch t, pscore(_pscore) outcome(y)
assert _rc==111
assert _pscore==p
capture cpsmatch t, pscore(p) neighbor(0)
assert _rc==111
capture cpsmatch t, pscore(p) kerneltype(unknown)
assert _rc==198
cpsmatch t, pscore(p) outcome(y) caliper(-1)
assert missing(r(att))
assert _nn==0 if t==1
cpsmatch t, pscore(p) outcome(y) kernel bwidth(0)
assert missing(r(att))
replace t=2*t
cpsmatch t, pscore(p) outcome(y)
assert missing(r(att))
assert missing(_nn)
assert _weight==0
noi di "OCT22_REGRESSIONS_COMPLETE"
    global CTOOLS_CASE_COMPLETE "oct22"
    exit
}

if "`__case'" == "oct22_reference" {
* Additional parity against installed reference ado commands.
* Install ivreg2/ranktest, binscatter and psmatch2 before running.
version 16
clear all
set more off
adopath ++ build
local candidate : environment CTOOLS_TEST_PLUGIN_DIR
if `"`candidate'"' != "" adopath ++ `"`candidate'"'
local reference : environment CTOOLS_TEST_REFERENCE_DIR
if `"`reference'"' != "" adopath ++ `"`reference'"'
set obs 100
set seed 2345
gen double z=rnormal()
gen double z2=rnormal()
gen double x=z+.5*z2+rnormal()
gen double y=10+.2*x+rnormal()
gen cl=ceil(_n/5)
foreach options in "" "small" "robust" "robust small" "cluster(cl)" "cluster(cl) small" ///
    "dofminus(2)" "small dofminus(2) sdofminus(3)" "robust dofminus(2) sdofminus(3)" ///
    "small robust dofminus(2) sdofminus(3)" "gmm2s robust" "gmm2s robust small" ///
    "gmm2s robust dofminus(2) sdofminus(3)" "gmm2s robust small dofminus(2) sdofminus(3)" ///
    "noconstant small" "cluster(cl) small dofminus(2) sdofminus(3)" {
    quietly ivreg2 y (x=z z2), `options'
    matrix expected=e(V)
    scalar expected_F=e(F)
    scalar expected_rmse=e(rmse)
    scalar expected_r2a=e(r2_a)
    quietly civreghdfe y (x=z z2), `options'
    matrix actual=e(V)
    assert mreldif(actual,expected)<1e-8
    assert reldif(e(F),expected_F)<1e-8
    assert reldif(e(rmse),expected_rmse)<1e-8
    assert reldif(e(r2_a),expected_r2a)<1e-8
}
* Effective sample size and partialled controls must retain the same corrections.
gen fw=1+mod(_n,3)
gen double aw=.2+runiform()
foreach weight in "" "[fw=fw]" "[aw=aw]" {
    foreach options in "partial(z2)" "partial(z2) small" ///
        "partial(z2) robust dofminus(2) sdofminus(3)" ///
        "partial(z2) cluster(cl) small dofminus(2) sdofminus(3)" {
        quietly ivreg2 y z2 (x=z) `weight', `options'
        matrix expected=e(V)
        scalar expected_F=e(F)
        scalar expected_rmse=e(rmse)
        quietly civreghdfe y z2 (x=z) `weight', `options'
        matrix actual=e(V)
        * Compare the common slope: ivreg2 partial() omits the constant from e(V).
        matrix actual=actual[1,1]
        assert mreldif(actual,expected)<1e-8
        assert reldif(e(F),expected_F)<1e-8
        assert reldif(e(rmse),expected_rmse)<1e-8
    }
}
* Covariance rank deficiency does not change ivreg2's slope-count divisor.
gen two_clusters=mod(_n,2)
gen double large_scale=1e6*rnormal()
gen double small_scale=1e-6*rnormal()
foreach rhs in "z2" "large_scale small_scale" {
    quietly ivreg2 y `rhs' (x=z), cluster(two_clusters) small
    scalar expected_F=e(F)
    scalar expected_df=e(df_m)
    quietly civreghdfe y `rhs' (x=z), cluster(two_clusters) small
    assert e(df_m)==expected_df
    assert reldif(e(F),expected_F)<1e-7
}
* Singleton fit requirements differ from retaining singleton scatter points.
clear
set obs 8
gen x=mod(_n,2)
gen y=_n
gen g=(_n==8)
foreach line in none lfit qfit {
    capture quietly binscatter y x, by(g) discrete linetype(`line') nodraw
    local expected_rc=_rc
    capture quietly cbinscatter y x, by(g) discrete linetype(`line') nograph
    assert _rc==`expected_rc'
}
* Reference treatment and option semantics, not assumptions about valid inputs.
foreach options in "neighbor(0)" "neighbor(-1)" "kernel bwidth(0)" "caliper(-1)" "kerneltype(unknown)" {
    clear
    set obs 10
    gen t=mod(_n,2)
    gen p=_n/11
    gen y=_n
    capture quietly psmatch2 t, pscore(p) outcome(y) `options'
    local expected_rc=_rc
    scalar expected_att=r(att)
    capture quietly cpsmatch t, pscore(p) outcome(y) `options'
    assert _rc==`expected_rc'
    if !`expected_rc' assert r(att)==expected_att
}
noi di "OCT22_REFERENCE_COMPLETE"
    global CTOOLS_CASE_COMPLETE "oct22_reference"
    exit
}

if "`__case'" == "oct27" {
/* Regression contracts for the HDFE wrapper fixes of 27 Oct 2026. */
do validation/support/validate_setup.do

capture noisily {
    sysuse auto, clear
    creghdfe price weight, absorb(rep78 foreign)
    assert e(N_hdfe) == 2 & e(N_hdfe_extended) == 2
    assert "`e(vcetype)'" == ""
    creghdfe price weight, absorb(rep78) vce(robust)
    assert "`e(vcetype)'" == "Robust" & "`e(vce)'" == "robust"
    creghdfe price weight, absorb(rep78) vce(cluster turn)
    assert "`e(vcetype)'" == "Robust" & "`e(vce)'" == "cluster"
    cppmlhdfe price weight, absorb(rep78)
    assert "`e(vcetype)'" == "Robust"
    cppmlhdfe price weight, absorb(rep78) vce(cluster turn)
    assert "`e(vcetype)'" == "Robust"
    * The first-stage regressions run by creghdfe carry the label as well
    civreghdfe price (weight = length), absorb(rep78) vce(robust) savefirst
    estimates restore _civreghdfe_weight
    assert "`e(vcetype)'" == "Robust"
    capture which reghdfe
    if !_rc {
        reghdfe price weight, absorb(rep78 foreign)
        assert e(N_hdfe_extended) == 2
    }
}
if _rc test_fail "N_hdfe_extended and Robust vcetype" "r(`=_rc')"
else test_pass "N_hdfe_extended and Robust vcetype"

capture noisily {
    sysuse auto, clear
    gen str5 srep = string(rep78) if !missing(rep78)
    gen str1 sfor = cond(foreign, "F", "D")
    gen int cell = 10 * rep78 + foreign
    creghdfe price weight, absorb(rep78 foreign) vce(cluster turn)
    matrix b0 = e(b)
    matrix V0 = e(V)
    local N0 = e(N)
    creghdfe price weight, absorb(srep sfor) vce(cluster turn)
    assert mreldif(b0, e(b)) < 1e-12 & mreldif(V0, e(V)) < 1e-12 & e(N) == `N0'
    assert "`e(absorb)'" == "srep sfor"
    creghdfe price weight, absorb(i.rep78 i.foreign) vce(cluster turn)
    assert mreldif(b0, e(b)) < 1e-12 & "`e(absorb)'" == "i.rep78 i.foreign"
    creghdfe price weight length, absorb(cell)
    matrix b1 = e(b)
    matrix V1 = e(V)
    creghdfe price weight length, absorb(srep#sfor)
    assert mreldif(b1, e(b)) < 1e-12 & mreldif(V1, e(V)) < 1e-12
    assert "`e(absorb)'" == "srep#sfor" & e(N_hdfe) == 1
    creghdfe price weight length, absorb(i.rep78#i.foreign)
    assert mreldif(b1, e(b)) < 1e-12 & "`e(absorb)'" == "i.rep78#i.foreign"
    capture which reghdfe
    if !_rc {
        reghdfe price weight length, absorb(rep78#foreign)
        assert mreldif(b1, e(b)) < 1e-10 & mreldif(V1, e(V)) < 1e-10
    }
    civreghdfe price (weight = length) turn, absorb(rep78 foreign) vce(cluster turn)
    matrix b2 = e(b)
    matrix V2 = e(V)
    civreghdfe price (weight = length) turn, absorb(srep sfor) vce(cluster turn)
    assert mreldif(b2, e(b)) < 1e-12 & mreldif(V2, e(V)) < 1e-12
    assert "`e(absorb)'" == "srep sfor"
    civreghdfe price (weight = length) turn, absorb(cell)
    matrix b3 = e(b)
    matrix V3 = e(V)
    civreghdfe price (weight = length) turn, absorb(srep#sfor)
    assert mreldif(b3, e(b)) < 1e-12 & mreldif(V3, e(V)) < 1e-12
    assert "`e(absorb)'" == "srep#sfor"
}
if _rc test_fail "string, i. and # absorbed effects" "r(`=_rc')"
else test_pass "string, i. and # absorbed effects"

capture noisily {
    sysuse auto, clear
    foreach bad in "c.length#i.rep78" "rep78##foreign" "L.rep78" "fe=rep78" {
        capture creghdfe price weight, absorb(`bad')
        assert _rc == 198
        capture civreghdfe price (weight = length), absorb(`bad')
        assert _rc == 198
    }
    * savefe is reghdfe's only absorb() suboption; civreghdfe takes none
    capture creghdfe price weight, absorb(rep78, keepsingletons)
    assert _rc == 198
    capture civreghdfe price (weight = length), absorb(rep78, savefe)
    assert _rc == 198
    creghdfe price weight, absorb(rep78, savefe)
    confirm numeric variable __hdfe1__
}
if _rc test_fail "unsupported absorb() terms and suboptions" "r(`=_rc')"
else test_pass "unsupported absorb() terms and suboptions"

capture noisily {
    sysuse auto, clear
    civreghdfe price (weight = length) if foreign == 0 [aw=mpg], absorb(rep78) vce(robust) small
    local cmdline `"`e(cmdline)'"'
    assert strpos(`"`cmdline'"', "civreghdfe price (weight = length)") == 1
    assert strpos(`"`cmdline'"', "if foreign == 0") > 0
    assert strpos(`"`cmdline'"', "[aw=mpg]") > 0
    assert strpos(`"`cmdline'"', "vce(robust) small") > 0
}
if _rc test_fail "civreghdfe e(cmdline) keeps if, weights and options" "r(`=_rc')"
else test_pass "civreghdfe e(cmdline) keeps if, weights and options"

capture noisily {
    sysuse auto, clear
    gen byte cl2 = mod(_n, 7)
    creghdfe price weight, absorb(rep78)
    foreach option in first rf savefirst saverf {
        capture civreghdfe price (weight = length), absorb(rep78) vce(cluster turn cl2) `option'
        assert _rc == 198
        assert "`e(cmd)'" == "creghdfe"
    }
    civreghdfe price (weight = length), absorb(rep78) vce(cluster turn cl2)
    assert "`e(vcetype)'" == "Two-way Cluster"
}
if _rc test_fail "two-way clustering refuses auxiliary regressions up front" "r(`=_rc')"
else test_pass "two-way clustering refuses auxiliary regressions up front"

capture noisily {
    sysuse auto, clear
    cppmlhdfe price weight, absorb(rep78 foreign) dof(pairwise clusters)
    matrix b0 = e(b)
    matrix V0 = e(V)
    local df_a = e(df_a)
    assert "`e(dofmethod)'" == "pairwise clusters"
    cppmlhdfe price weight, absorb(rep78 foreign) dof(pair clus)
    assert mreldif(b0, e(b)) == 0 & mreldif(V0, e(V)) == 0 & e(df_a) == `df_a'
    assert "`e(dofmethod)'" == "pairwise clusters"
    cppmlhdfe price weight, absorb(rep78 foreign) dof(firstp)
    assert "`e(dofmethod)'" == "firstpair"
    cppmlhdfe price weight, absorb(rep78 foreign)
    assert "`e(dofmethod)'" == "pairwise clusters continuous"
    foreach bad in "pair first" "all clusters" "junk" {
        capture cppmlhdfe price weight, absorb(rep78 foreign) dof(`bad')
        assert _rc == 198
    }
    creghdfe price weight, absorb(rep78 foreign) dof(pair clus)
}
if _rc test_fail "dofadjustments() abbreviations" "r(`=_rc')"
else test_pass "dofadjustments() abbreviations"

capture noisily {
    sysuse auto, clear
    civreghdfe price (weight = length), absorb(rep78) noreturn
    mata: st_local("left", invtokens((st_dir("global", "numscalar", "__civreghdfe_*") \ ///
        st_dir("global", "matrix", "__civreghdfe_*"))'))
    assert "`left'" == ""
}
if _rc test_fail "civreghdfe noreturn drops its plugin scalars" "r(`=_rc')"
else test_pass "civreghdfe noreturn drops its plugin scalars"

capture noisily {
    sysuse auto, clear
    gen byte cl2 = mod(_n, 7)
    gen str3 smk = substr(make, 1, 3)
    gen int cell = 100 * turn + foreign
    * String, interaction and i. cluster dimensions, as reghdfe names them
    creghdfe price weight, absorb(rep78) vce(cluster cell)
    matrix b0 = e(b)
    matrix V0 = e(V)
    local G0 = e(N_clust)
    creghdfe price weight, absorb(rep78) vce(cluster turn#foreign)
    assert mreldif(b0, e(b)) < 1e-12 & mreldif(V0, e(V)) < 1e-12 & e(N_clust) == `G0'
    assert "`e(clustvar)'" == "turn#foreign" & "`e(clustvar1)'" == "turn#foreign"
    creghdfe price weight, absorb(rep78) vce(cluster i.turn)
    assert "`e(clustvar)'" == "turn" & e(N_clustervars) == 1 & e(N_clust1) == e(N_clust)
    creghdfe price weight, absorb(rep78) vce(cluster smk)
    assert "`e(clustvar)'" == "smk" & e(N_clust) == 22
    matrix bs = e(b)
    matrix Vs = e(V)
    civreghdfe price (weight = length), absorb(rep78) vce(cluster cell)
    matrix b1 = e(b)
    matrix V1 = e(V)
    civreghdfe price (weight = length), absorb(rep78) vce(cluster turn#foreign)
    assert mreldif(b1, e(b)) < 1e-12 & mreldif(V1, e(V)) < 1e-12
    civreghdfe price (weight = length), absorb(rep78) vce(cluster turn cl2)
    assert "`e(clustvar)'" == "turn cl2" & "`e(clustvar1)'" == "turn" & "`e(clustvar2)'" == "cl2"
    cppmlhdfe price weight, absorb(rep78) vce(cluster i.turn foreign#rep78)
    assert "`e(clustvar)'" == "turn foreign#rep78" & e(N_clustervars) == 2
    assert e(N_clust) == min(e(N_clust1), e(N_clust2))
    foreach bad in "L.turn" "c.turn" {
        capture creghdfe price weight, absorb(rep78) vce(cluster `bad')
        assert _rc == 198
    }
    * pweights need robust or clustered standard errors, as in reghdfe
    capture creghdfe price weight [pw=mpg], absorb(rep78) vce(ols)
    assert _rc == 9
    creghdfe price weight [pw=mpg], absorb(rep78)
    assert "`e(vce)'" == "robust"
    capture which reghdfe
    if !_rc {
        reghdfe price weight, absorb(rep78) vce(cluster smk)
        assert mreldif(e(b), bs) < 1e-10 & mreldif(e(V), Vs) < 1e-10 & e(N_clust) == 22
    }
}
if _rc test_fail "shared vce() parsing" "r(`=_rc')"
else test_pass "shared vce() parsing"

capture noisily {
    sysuse auto, clear
    creghdfe price weight length, absorb(i.rep78 foreign#turn) dof(firstpair clusters) vce(cluster turn)
    assert "`e(absvars)'" == "i.rep78 foreign#turn"
    assert "`e(extended_absvars)'" == "rep78 foreign#turn"
    assert "`e(title)'" == "HDFE Linear regression" & "`e(title2)'" == "Absorbing 2 HDFE groups"
    assert "`e(dofmethod)'" == "firstpair clusters"
    assert e(N_hdfe_extended) == 2 & e(drop_singletons) == 1
    matrix t = e(dof_table)
    assert rowsof(t) == 2 & colsof(t) == 5
    assert e(df_a_initial) == t[1,1] + t[2,1] & e(df_a_redundant) == t[1,2] + t[2,2]
    assert t[2,5] == 1 & missing(t[1,4])
    local c_initial = e(df_a_initial)
    local c_redundant = e(df_a_redundant)
    local c_df_a = e(df_a)
    creghdfe price weight
    assert "`e(absvars)'" == "_cons" & e(df_a_initial) == 1 & e(df_a_redundant) == 0
    cppmlhdfe price weight, absorb(rep78 foreign)
    matrix t = e(dof_table)
    local columns : colnames t
    assert "`columns'" == "Categories Redundant Num_Coefs Inexact Nested"
    assert e(df_a_initial) == t[1,1] + t[2,1] & e(df_a_redundant) == t[1,2] + t[2,2]
    capture which reghdfe
    if !_rc {
        reghdfe price weight length, absorb(i.rep78 foreign#turn) dof(firstpair clusters) vce(cluster turn)
        assert e(df_a_initial) == `c_initial' & e(df_a_redundant) == `c_redundant' & e(df_a) == `c_df_a'
        * Cluster variables are marked out before the factor variables are
        * expanded, so both commands name the same coefficients
        gen cl3 = turn if foreign == 0
        reghdfe price weight i.foreign, absorb(rep78) vce(cluster cl3)
        local names_r : colfullnames e(b)
        creghdfe price weight i.foreign, absorb(rep78) vce(cluster cl3)
        local names_c : colfullnames e(b)
        assert "`names_r'" == "`names_c'"
    }
}
if _rc test_fail "shared reghdfe e() results" "r(`=_rc')"
else test_pass "shared reghdfe e() results"

print_summary oct27
if $TESTS_FAILED > 0 exit 9
global CTOOLS_COMPONENT_COMPLETE "oct27"
    global CTOOLS_CASE_COMPLETE "oct27"
    exit
}

if "`__case'" == "sep22" {
/* Regression contracts for SEP22_ASTRA_AUDIT.md. No network data required. */
do validation/support/validate_setup.do

capture noisily {
    tempfile using_data
    foreach relation in 1:1 m:1 1:m {
        foreach side in master using {
            if ("`relation'" == "m:1" & "`side'" == "master") | ("`relation'" == "1:m" & "`side'" == "using") continue
            clear
            set obs 3
            gen double id = _n
            if "`side'" == "using" replace id = 99 in 2/3
            gen str8 us = "using"
            save `using_data', replace
            clear
            set obs 3
            gen double id = _n
            if "`side'" == "master" replace id = 99 in 2/3
            gen long row = _n
            quietly datasignature
            local before "`r(datasignature)'"
            capture cmerge `relation' id using `using_data'
            assert _rc == 459
            quietly datasignature
            assert "`r(datasignature)'" == "`before'"
        }
    }
    * Missing string keys obey the same cardinality contract.
    clear
    set obs 2
    gen str5 key = ""
    gen row = _n
    save `using_data', replace
    capture cmerge 1:1 key using `using_data'
    assert _rc == 459
    assert _N == 2 & row == _n
}
if _rc test_fail "S01 merge cardinality and rollback" "r(`=_rc')"
else test_pass "S01 merge cardinality and rollback"

capture noisily {
    sysuse auto, clear
    gen str8 group = cond(foreign, "foreign", "domestic")
    clonevar price_before = price
    cwinsor price, suffix(_w) by(foreign group)
    assert price == price_before
    assert !missing(price_w)
    drop price_w
    order price, last
    cwinsor price if mpg > 20, suffix(_w)
    assert missing(price_w) if mpg <= 20
    cwinsor price, replace
    assert !missing(price)
    * A later invalid target cannot leave an earlier generated variable behind.
    drop price_w
    gen mpg_w = 1
    capture cwinsor price mpg, suffix(_w)
    assert _rc == 110
    capture confirm variable price_w
    assert _rc == 111
    assert mpg_w == 1
}
if _rc test_fail "S02 winsor destinations and rollback" "r(`=_rc')"
else test_pass "S02 winsor destinations and rollback"

capture noisily {
    clear
    set obs 400
    set seed 45678
    gen g = ceil(_n/20)
    gen double cl = g + cond(mod(_n,2), .1, .2)
    egen cid = group(cl)
    gen double x = rnormal()
    gen double y = .4*x + g/10 + rnormal()
    gen double yp = rpoisson(exp(.3*x+g/100))
    foreach cmd in creghdfe cppmlhdfe {
        local dep = cond("`cmd'" == "creghdfe", "y", "yp")
        `cmd' `dep' x, absorb(g) vce(cluster cl)
        assert e(N_clust) == 40
        matrix V = e(V)
        `cmd' `dep' x, absorb(g) vce(cluster cid)
        assert e(N_clust) == 40
        mata: assert(mreldif(st_matrix("V"),st_matrix("e(V)")) < 1e-7)
    }
}
if _rc test_fail "S03 cluster partition invariance" "r(`=_rc')"
else test_pass "S03 cluster partition invariance"

capture noisily {
    clear
    set obs 400
    set seed 45679
    gen g = ceil(_n/20)
    gen double frac = g/100
    gen double large = 1e12+g
    gen double neg = -g/100
    gen str8 strcl = string(g)
    gen h = mod(_n,17)
    gen double hfrac = h/100
    gen double z = rnormal()
    gen double x = z + rnormal()
    gen double y = .4*x + g/10 + rnormal()
    civreghdfe y (x=z), absorb(g) vce(cluster g)
    matrix V = e(V)
    foreach cl in frac large neg strcl {
        civreghdfe y (x=z), absorb(g) vce(cluster `cl')
        mata: assert(mreldif(st_matrix("V"),st_matrix("e(V)")) < 1e-7)
    }
    civreghdfe y (x=z), absorb(g) vce(cluster g h)
    matrix V = e(V)
    civreghdfe y (x=z), absorb(g) vce(cluster frac hfrac)
    mata: assert(mreldif(st_matrix("V"),st_matrix("e(V)")) < 1e-7)
    gen byte one_cluster = 1
    capture civreghdfe y (x=z), absorb(g) vce(cluster one_cluster)
    assert _rc == 459
    capture civreghdfe y (x=z), absorb(g) vce(cluster g one_cluster)
    assert _rc == 459
}
if _rc test_fail "S04 IV cluster labels" "r(`=_rc')"
else test_pass "S04 IV cluster labels"

capture noisily {
    clear
    set obs 402
    set seed 45680
    gen g = ceil(_n/20)
    replace g = 99 in 401
    replace g = 98 in 402
    gen double x = rnormal()
    gen double y = .4*x + g/10 + rnormal()
    gen double yp = rpoisson(exp(.3*x))
    creghdfe y x if _n <= 401, absorb(g) residuals(res)
    assert e(N) == 400
    count if e(sample)
    assert r(N) == 400
    assert e(sample) == (_n <= 400)
    assert !missing(res) == e(sample)
    estimates store ols_sample
    cppmlhdfe yp x if _n <= 401, absorb(g)
    assert e(N) == 400
    assert e(sample) == (_n <= 400)
    estimates restore ols_sample
    assert e(sample) == (_n <= 400)
}
if _rc test_fail "S05 final sample and residuals" "r(`=_rc')"
else test_pass "S05 final sample and residuals"

capture noisily {
    foreach zero in 1 10 20 {
        clear
        set obs 402
        set seed 45681
        gen g = ceil(_n/20)
        replace g = 99 in 401
        gen double cl = mod(_n,37)/100
        gen double x = rnormal()
        gen double off = rnormal()/10
        gen double w = 1+runiform()
        gen double y = rpoisson(exp(.3*x+off))
        replace y = 0 if g == `zero'
        cppmlhdfe y x if _n <= 401 [aw=w], absorb(g) vce(cluster cl) offset(off)
        assert e(N) == 380
        assert e(sample) == (_n <= 400 & g != `zero')
        matrix B = e(b)
        matrix V = e(V)
        local nc = e(N_clust)
        drop if _n > 400 | g == `zero'
        cppmlhdfe y x [aw=w], absorb(g) vce(cluster cl) offset(off)
        assert e(N_clust) == `nc'
        mata: assert(mreldif(st_matrix("B"),st_matrix("e(b)")) < 1e-7)
        mata: assert(mreldif(st_matrix("V"),st_matrix("e(V)")) < 1e-7)
    }
}
if _rc test_fail "S06 PPML separation and row alignment" "r(`=_rc')"
else test_pass "S06 PPML separation and row alignment"

capture noisily {
    clear
    set obs 200
    set seed 922
    gen g = ceil(_n/20)
    gen double x = rnormal()
    gen y = (_n < 110)
    gen byte one = 1
    capture cppmlhdfe y x, absorb(g) vce(cluster one)
    assert _rc == 459
    * A deliberately large guard threshold must fail, not count retained rows
    * as repeatedly dropped. A subsequent ordinary fit must remain usable.
    capture cppmlhdfe y x, absorb(g) septolerance(10)
    assert _rc == 0
    cppmlhdfe y x, absorb(g)
    assert e(num_singletons) == 80
    assert e(num_separated) == 0
    assert e(sample) == (g <= 6)
}
if _rc test_fail "S05/S06 PPML guards and recovery" "r(`=_rc')"
else test_pass "S05/S06 PPML guards and recovery"

capture noisily {
    clear
    set obs 8
    gen double x = _n
    replace x = .a in 7
    replace x = .z in 8
    mata: st_vlmodify("lab", (1\2\3\4\5\6\.a\.z), ("plain"\(char(96)+"literal"+char(39))\("quote "+char(34)+" text")\("tabs"+char(9)+"line"+char(10)+char(13))\"café 文"\(3000*"x")\"refused"\"unknown"))
    label values x lab
    decode x, gen(native)
    cdecode x, gen(actual)
    assert r(N_vars) == 1
    assert actual == native
    assert strlen(actual[6]) == 3000
    cdecode x, gen(short) maxlength(4)
    decode x, gen(native_short) maxlength(4)
    assert short == native_short
    cdecode x, replace
    assert x == native
    * All variables are staged before replacing any input.
    gen y = 1
    label values y lab
    gen unlabeled = 2
    capture cdecode y unlabeled, replace
    assert _rc == 182
    confirm numeric variable y
}
if _rc test_fail "S07 literal labels and long strings" "r(`=_rc')"
else test_pass "S07 literal labels and long strings"

capture noisily {
    clear
    set obs 40
    gen cl = ceil(_n/10)
    gen str3 st = cond(cl<=2,"a","b")
    gen freq = 0
    cbsample 1, cluster(cl) weight(freq)
    summarize freq, meanonly
    assert r(sum) == 10
    cbsample 2, strata(cl) weight(freq)
    summarize freq, meanonly
    assert r(sum) == 8
    cbsample, strata(cl) weight(freq)
    summarize freq, meanonly
    assert r(sum) == 40
    cbsample 1, strata(st) cluster(cl) weight(freq)
    bysort st: egen drawn = total(freq)
    assert drawn == 10
    capture cbsample 3, strata(st) cluster(cl) weight(freq)
    assert _rc == 498
    assert _N == 40
    drop drawn
    * Repeated cluster labels in different strata are separate sampling units.
    replace cl = mod(cl,2)
    cbsample 1, strata(st) cluster(cl) weight(freq)
    bysort st: egen drawn = total(freq)
    assert drawn == 10
    cbsample 2, strata(st)
    assert _N == 4
}
if _rc test_fail "S08 bootstrap draw counts" "r(`=_rc')"
else test_pass "S08 bootstrap draw counts"

capture noisily {
    clear
    set obs 400
    set seed 45682
    gen g = ceil(_n/20)
    gen double x = rnormal()
    gen double z = x+rnormal()
    gen double y = .4*x + rnormal()
    gen double yp = rpoisson(exp(.3*x))
    gen double w = 1+runiform()
    gen double twice = 2*w
    foreach cmd in creghdfe cppmlhdfe civreghdfe {
        local dep = cond("`cmd'" == "cppmlhdfe", "yp", "y")
        local rhs = cond("`cmd'" == "civreghdfe", "(x=z)", "x")
        `cmd' `dep' `rhs' [aw=2*w], absorb(g)
        matrix B = e(b)
        matrix V = e(V)
        `cmd' `dep' `rhs' [aw=twice], absorb(g)
        mata: assert(mreldif(st_matrix("B"),st_matrix("e(b)")) < 1e-7)
        mata: assert(mreldif(st_matrix("V"),st_matrix("e(V)")) < 1e-7)
    }
    cppmlhdfe yp x, absorb(g)
    estimates store ppml_prediction
    capture predict double mu
    assert _rc == 198
    capture confirm variable mu
    assert _rc == 111
    predict double xb, xb
    assert reldif(xb, x*_b[x]+_b[_cons]) < 1e-7
    estimates restore ppml_prediction
    capture predict double mu
    assert _rc == 198
}
if _rc test_fail "S09/S10 predictions and weight expressions" "r(`=_rc')"
else test_pass "S09/S10 predictions and weight expressions"

capture noisily {
    clear
    input byte treat double(ps y)
    0 .5 1
    0 .5 3
    1 .5 7
    end
    cpsmatch treat, pscore(ps) outcome(y) ties
    assert _nn == 2 if treat
    assert missing(_nn) if !treat
    cpsmatch treat, pscore(ps) outcome(y) neighbor(5)
    assert _nn == 2 if treat
    cpsmatch treat, pscore(ps) outcome(y) radius caliper(.1)
    assert _nn == 2 if treat
    cpsmatch treat, pscore(ps) outcome(y) kernel bwidth(.1)
    assert _nn == 2 if treat
}
if _rc test_fail "S11 actual matching neighbor counts" "r(`=_rc')"
else test_pass "S11 actual matching neighbor counts"

capture noisily {
    tempfile empty full
    clear
    set obs 2
    gen id = _n
    gen str10 text = "example"
    save `full', replace
    drop in 1/2
    save `empty', replace
    foreach nogenerate in "" "nogenerate" {
        use `full', clear
        cmerge 1:1 id using `empty', keep(match) `nogenerate'
        assert _N == 0
        use `empty', clear
        cmerge 1:1 id using `full', keep(master) `nogenerate'
        assert _N == 0
        use `empty', clear
        cmerge 1:1 id using `full', keep(using) `nogenerate'
        assert _N == 2
        local typ : type text
        assert "`typ'" == "str10"
    }
    sysuse auto, clear
    capture cbinscatter price mpg, genxq(q) nograph
    assert _rc == 198
    capture confirm variable q
    assert _rc == 111
}
if _rc test_fail "S12/S15 empty-side filtering and unsupported option" "r(`=_rc')"
else test_pass "S12/S15 empty-side filtering and unsupported option"

print_summary sep22
if $TESTS_FAILED > 0 exit 9

* Reached only after the complete component script.
global CTOOLS_COMPONENT_COMPLETE "sep22"
    global CTOOLS_CASE_COMPLETE "sep22"
    exit
}

if "`__case'" == "sep24" {
/* Regression contracts for the 21 SEP24_ASTRA_AUDIT findings. */
do validation/support/validate_setup.do

capture noisily {
    tempfile using_data
    clear
    set obs 2
    gen long id = 1000 + _n
    gen double value = 1000.125
    gen str20 text = "abcdefghijklmnopqrst"
    save `using_data', replace
    foreach option in "" "update" "update replace" {
        clear
        set obs 1
        gen byte id = 1
        gen byte value = .
        gen str1 text = ""
        cmerge 1:1 id using `using_data', `option'
        assert id == 1000 + _n - 1 if _n > 1
        assert value == 1000.125 if _merge == 2
        assert text == "abcdefghijklmnopqrst" if _merge == 2
        local vt : type value
        local st : type text
        assert "`vt'" == "double" & "`st'" == "str20"
    }
    clear
    set obs 1
    gen str4 id = "1001"
    gen row = _n
    quietly datasignature
    local before "`r(datasignature)'"
    foreach option in "" "force" {
        capture cmerge 1:1 id using `using_data', `option'
        assert _rc == 106
        quietly datasignature
        assert "`before'" == "`r(datasignature)'"
    }
}
if _rc test_fail "A01/A08 merge promotion and type safety" "r(`=_rc')"
else test_pass "A01/A08 merge promotion and type safety"

capture noisily {
    sysuse auto, clear
    quietly datasignature
    local before "`r(datasignature)'"
    capture creghdfe price mpg, residuals(mpg)
    assert _rc == 110
    capture crangestat (mean) mpg=mpg, interval(weight -100 100)
    assert _rc == 110
    capture crangestat (mean) out=mpg (sum) out=price, interval(weight -100 100)
    assert _rc == 198
    quietly datasignature
    assert "`before'" == "`r(datasignature)'"
    capture confirm variable out
    assert _rc == 111
}
if _rc test_fail "A02 output aliases and duplicates preserve data" "r(`=_rc')"
else test_pass "A02 output aliases and duplicates preserve data"

capture noisily {
    clear
    set obs 1000
    set seed 296519
    gen double x = rnormal()
    gen double z = x+rnormal()
    gen double z2 = 2*z
    gen double tiny_z = 1e-9*z
    gen double y = .3*x+4*z+rnormal()
    gen g = mod(_n,10)
    foreach absorb in "" "absorb(g)" {
        cbinscatter y x, controls(z) method(binsreg) nograph `absorb'
        matrix one = e(bindata)
        cbinscatter y x, controls(tiny_z) method(binsreg) nograph `absorb'
        matrix scaled = e(bindata)
        mata: assert(max(abs(st_matrix("one")-st_matrix("scaled"))) < 1e-7)
        cbinscatter y x, controls(z z2) method(binsreg) nograph `absorb'
        matrix dup = e(bindata)
        mata: assert(max(abs(st_matrix("one")-st_matrix("dup"))) < 1e-7)
    }
    gen double w = 1+abs(z)
    cbinscatter y x [aw=w], nograph
    matrix one = e(bindata)
    cbinscatter y x [aw=1+abs(z)], nograph
    matrix dup = e(bindata)
    mata: assert(max(abs(st_matrix("one")-st_matrix("dup"))) < 1e-10)
    gen byte category = mod(_n,4)
    quietly tab category, gen(dummy)
    cbinscatter y x, controls(dummy2 dummy3 dummy4) nograph
    matrix one = e(bindata)
    cbinscatter y x, controls(i.category) nograph
    matrix dup = e(bindata)
    mata: assert(max(abs(st_matrix("one")-st_matrix("dup"))) < 1e-7)
}
if _rc test_fail "A03/A16 redundant and factor controls, weight expressions" "r(`=_rc')"
else test_pass "A03/A16 redundant and factor controls, weight expressions"

capture noisily {
    clear
    set obs 1000
    set seed 742913
    gen g = ceil(_n/10)
    gen h = ceil(runiform()*50)
    gen double z = rnormal()
    gen double x = z+rnormal()+g/10+h/5
    gen double y = .7*x+g/5-h/7+rnormal()
    capture creghdfe y x, absorb(g h) iterate(1) tolerance(1e-14) residuals(res)
    assert _rc == 430
    capture confirm variable res
    assert _rc == 111
    capture civreghdfe y (x=z), absorb(g h) maxiter(1) tolerance(1e-14)
    assert _rc == 430
    foreach option in "iterate(0)" "tolerance(0)" {
        capture creghdfe y x, absorb(g h) `option'
        assert _rc == 198
    }
    creghdfe y x, absorb(g h) tolerance(1e-10)
    assert !missing(_b[x])
    replace y = rpoisson(exp(.3+.05*x+g/200+h/100))
    capture cppmlhdfe y x, absorb(g h) irlsmaxiter(1)
    assert _rc == 430
}
if _rc test_fail "A04 exhausted solvers return errors" "r(`=_rc')"
else test_pass "A04 exhausted solvers return errors"

capture noisily {
    foreach n in 999 1000 1001 2001 {
        clear
        set obs `n'
        gen double key = _n
        replace key = .a in `n'
        gen double v = 1
        replace v = 1001 in `n'
        foreach excl in "" "excludeself" {
            crangestat (mean) actual=v (sum) total=v (count) count=v, interval(key . .) `excl'
            assert actual == 1 if !missing(key)
            assert total == `n'-1-("`excl'"!="") if !missing(key)
            assert count == total if !missing(key)
            drop actual total count
        }
    }
}
if _rc test_fail "A07 missing range keys across thresholds" "r(`=_rc')"
else test_pass "A07 missing range keys across thresholds"

capture noisily {
    foreach n in 31 32 100 30000 {
        clear
        set obs `n'
        gen str8 key = cond(mod(_n,3)==0,"",cond(mod(_n,2),"az","ba"))
        gen long row = _n
        csort key, algorithm(merge) nosortedby threads(2)
        assert key >= key[_n-1] if _n > 1
        assert row > row[_n-1] if _n > 1 & key == key[_n-1]
    }
}
if _rc test_fail "A12 merge lexicographic order and stable ties" "r(`=_rc')"
else test_pass "A12 merge lexicographic order and stable ties"

capture noisily {
    clear
    set obs 600
    set seed 92815
    gen g = ceil(_n/10)
    gen double x = rnormal()
    gen double z = rnormal()
    gen double zero = 0
    gen double y = rpoisson(exp(.3+.4*x+.2*z+g/200))
    gen byte small_cluster = 1+mod(g,2)
    foreach vc in "vce(robust)" "vce(cluster g)" "vce(cluster small_cluster)" {
    foreach vars in "zero x z" "x zero z" "x z zero" {
        ppmlhdfe y `vars', absorb(g) `vc'
        scalar reference_chi2 = e(chi2)
        scalar reference_df = e(df_m)
        cppmlhdfe y `vars', absorb(g) `vc'
        assert e(df_m) == reference_df
        if missing(reference_chi2) assert missing(e(chi2))
        else assert reldif(e(chi2),reference_chi2) < 1e-6
    }
    }
}
if _rc test_fail "A14 PPML joint tests select retained regressors" "r(`=_rc')"
else test_pass "A14 PPML joint tests select retained regressors"

capture noisily {
    clear
    set obs 1
    gen str1 a = "x"
    gen strL b = "z"
    forvalues i=1/12 {
        quietly replace b=b+b
    }
    quietly datasignature
    local before "`r(datasignature)'"
    capture cencode a b, replace
    assert _rc != 0
    quietly datasignature
    assert "`before'" == "`r(datasignature)'"
    label define first 1 "old"
    capture cencode a b, gen(first second)
    assert _rc != 0
    local old : label first 1
    local added : label first 2, strict
    assert "`old'" == "old" & "`added'" == ""
    gen byte target=9
    drop in 1
    capture cencode a, gen(target)
    assert _rc == 110
    local vt : type target
    assert "`vt'" == "byte"
    capture cencode a b, gen(new new)
    assert _rc == 198
    capture cencode target, gen(new)
    assert _rc == 107
    cencode a, gen(new)
    confirm numeric variable new
}
if _rc test_fail "A15 encoding failure and empty data contracts" "r(`=_rc')"
else test_pass "A15 encoding failure and empty data contracts"

capture noisily {
    clear
    set obs 3
    gen double day = _n-1
    format day %td
    gen double timestamp = clock("24sep2026 12:34:56.789","DMYhms")+_n
    format timestamp %tc
    gen str10 note = "hello " + string(_n)
    foreach ext in csv xlsx {
        local path "temp/space threads(2) name.`ext'"
        if "`ext'" == "csv" cexport delimited using "`path'", replace
        else cexport excel using "`path'", firstrow(variables) sheet("A sheet") replace
        confirm file "`path'"
        if "`ext'" == "csv" {
            capture cexport delimited using "`path'"
        }
        else {
            capture cexport excel using "`path'"
        }
        assert _rc == 602
    }
    cimport excel using "temp/space threads(2) name.xlsx", firstrow clear
    assert day == _n-1
    assert note == "hello " + string(_n)
    assert abs(timestamp - (clock("24sep2026 12:34:56.789","DMYhms")+_n)) < .1
    cimport excel using "validation/fixtures/sep24_dates1900.xlsx", firstrow clear
    assert day == td(28feb1900) in 1
    * Native Excel import maps fictitious 29feb1900 (serial 60) to 28feb1900.
    assert day == td(28feb1900) in 2
    assert day == td(01mar1900) in 3
    assert day == 0 in 4
    cimport excel using "validation/fixtures/sep24_dates1904.xlsx", firstrow clear
    assert day == td(01jan1904) in 1
    assert day == 0 in 2
}
if _rc test_fail "A05/A06 dates, datetimes, literal export paths" "r(`=_rc')"
else test_pass "A05/A06 dates, datetimes, literal export paths"
capture erase "temp/space threads(2) name.csv"
capture erase "temp/space threads(2) name.xlsx"

capture noisily {
    clear
    set obs 1
    forvalues i=1/1100 {
        gen double long_variable_name_for_test_`i' = `i'+.125
    }
    cexport delimited using "temp/sep24_wide.csv", replace
    cexport excel using "temp/sep24_wide.xlsx", firstrow(variables) replace
    import delimited using "temp/sep24_wide.csv", clear varnames(1)
    assert c(k) == 1100
    forvalues i=1/1100 {
        assert long_variable_name_for_test_`i' == `i'+.125
    }
    import excel using "temp/sep24_wide.xlsx", clear firstrow
    assert c(k) == 1100
    forvalues i=1/1100 {
        assert long_variable_name_for_test_`i' == `i'+.125
    }
}
if _rc test_fail "A13 complete wide export metadata" "r(`=_rc')"
else test_pass "A13 complete wide export metadata"
capture erase "temp/sep24_wide.csv"
capture erase "temp/sep24_wide.xlsx"

capture noisily {
    _ctools_load
    assert "`r(plugin_version)'" != ""
    local revision "`r(plugin_revision)'"
    _ctools_load
    assert "`r(plugin_revision)'" == "`revision'"
}
if _rc test_fail "A20 repeated shared plugin loading" "r(`=_rc')"
else test_pass "A20 repeated shared plugin loading"

print_summary sep24
if $TESTS_FAILED > 0 exit 9
global CTOOLS_COMPONENT_COMPLETE "sep24"
    global CTOOLS_CASE_COMPLETE "sep24"
    exit
}

di as error "Unknown regressions case: `__case'"
exit 198
