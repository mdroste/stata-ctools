* s2_parquet: named cases share this file, each do-call has its own local scope.
gettoken __case 0 : 0
macro shift
if "`__case'" == "" local __case "s2_parquet"

if "`__case'" == "s2_parquet" {
do validation/support/validate_setup.do
set seed 20261104
* Fixtures: validation/support/parquet_fixtures.py writes them.
local fx "validation/fixtures/parquet"
* The reference, StataNow 19.5's import parquet, must be available.
capture findfile import_parquet.ado
if _rc {
    test_fail "IO-01 reference" "import parquet (StataNow 19.5) not found"
    print_summary s2_parquet
    exit 9
}

capture program drop _pq_run
program define _pq_run, rclass
    * Runs one command with its output logged; returns rc, r() and the log
    args command log
    label drop _all
    quietly log using `"`log'"', text replace name(_pq_log)
    capture noisily `command'
    local rc = _rc
    local N = r(N)
    local k = r(k)
    local Nr = r(N_rows)
    local nc = r(n_columns)
    quietly log close _pq_log
    return scalar rc = `rc'
    return scalar N = `N'
    return scalar k = `k'
    return scalar N_rows = `Nr'
    return scalar n_columns = `nc'
end

capture program drop pq_compare
program define pq_compare
    * import parquet and cimport parquet on the same arguments: return code,
    * stored results, variables (names, types, formats, labels, value
    * labels), values (cf), and the messages shown
    args label arguments setup
    tempfile ref reflog newlog
    if `"`setup'"' != "" quietly `setup'
    else clear
    _pq_run `"import parquet `arguments'"' `"`reflog'"'
    local rc1 = r(rc)
    local results1 = "`r(N)' `r(k)' `r(N_rows)' `r(n_columns)'"
    mata: st_global("PQ_SIG_ref", _pq_signature())
    if !`rc1' & c(k) quietly save `"`ref'"', replace
    if `"`setup'"' != "" quietly `setup'
    else clear
    _pq_run `"cimport parquet `arguments'"' `"`newlog'"'
    local rc2 = r(rc)
    local results2 = "`r(N)' `r(k)' `r(N_rows)' `r(n_columns)'"
    mata: st_global("PQ_SIG_new", _pq_signature())
    local why ""
    if `rc1' != `rc2' local why "rc `rc1' vs `rc2'"
    else if "`results1'" != "`results2'" local why "r() `results1' vs `results2'"
    else {
        mata: st_local("same", strofreal(st_global("PQ_SIG_ref") == st_global("PQ_SIG_new")))
        if !`same' local why "variables differ"
        else if !`rc1' & c(k) {
            capture cf _all using `"`ref'"', all
            if _rc local why "values differ"
        }
    }
    if "`why'" == "" {
        mata: st_local("same", strofreal(_pq_same_output(st_local("reflog"), st_local("newlog"))))
        if !`same' local why "output differs"
    }
    if "`why'" == "" test_pass `"`label'"'
    else test_fail `"`label'"' `"`why'"'
    macro drop PQ_SIG_*
end

mata:
/* Variables and value labels, read in Mata so that file-made names and
   labels never pass through macro expansion */
string scalar _pq_signature()
{
    real scalar j
    string scalar s, lab
    real colvector values
    string colvector text
    s = ""
    for (j = 1; j <= st_nvar(); j++) {
        s = s + "|" + st_varname(j) + ":" + st_vartype(j) + ":" + st_varformat(j) + ":" +
            st_varlabel(j) + ":" + st_varvaluelabel(j)
        lab = st_varvaluelabel(j)
        if (lab != "") {
            st_vlload(lab, values, text)
            if (rows(values)) s = s + "{" + invtokens(strofreal(values', "%21.0g") :+ "=" :+ text', char(1)) + "}"
        }
    }
    return(s)
}

/* The logged output of each command, without the Mata traceback lines that
   name Stata's internal functions: a traceback line's error message is
   kept (Stata cuts it at 80 columns, so it need only begin the other's) */
string colvector _pq_lines(string scalar path, real colvector traced)
{
    real scalar fh
    string scalar line
    string colvector out
    out = J(0, 1, "")
    traced = J(0, 1, .)
    fh = fopen(path, "r")
    while ((line = fget(fh)) != J(0, 0, "")) {
        line = strrtrim(line)
        if (strtrim(line) == "") continue
        if (regexm(line, "^ *([A-Za-z_0-9]+\(\)|<istmt>): +([0-9]+|-) +(.*)$")) {
            if (regexs(2) == "-" | regexs(3) == "Stata returned error") continue
            out = out \ regexs(3)
            traced = traced \ 1
            continue
        }
        out = out \ line
        traced = traced \ 0
    }
    fclose(fh)
    return(out)
}

real scalar _pq_same_output(string scalar reflog, string scalar newlog)
{
    real colvector t1, t2
    string colvector a, b
    real scalar i
    a = _pq_lines(reflog, t1)
    b = _pq_lines(newlog, t2)
    if (rows(a) != rows(b)) return(0)
    for (i = 1; i <= rows(a); i++) {
        if (a[i] == b[i]) continue
        if (t1[i] & strpos(b[i], a[i]) == 1) continue
        return(0)
    }
    return(1)
}
end

* ---------------------------------------------------------------------------
* Every fixture: import and describe
* ---------------------------------------------------------------------------
local files : dir "`fx'" files "*.parquet"
foreach f of local files {
    pq_compare "IO-01 fixture `f' import" `"using "`fx'/`f'", clear"'
    pq_compare "IO-01 fixture `f' describe" `"using "`fx'/`f'", describe"'
}

* ---------------------------------------------------------------------------
* Options
* ---------------------------------------------------------------------------
local j 0
foreach r in "2:4" "4" ":3" "-3" "-3:-1" "f:l" "F:L" "3:" "5:4" "0" "11" "0:3" "2:11" "1.5" ///
    "a" "2:b" "  2 : 4 " "-10" "-11" "f" "l" "3:l" ":l" ":-2" "2:0" "1:2.5" "1.0:2" "2:3.0" ///
    "-0" "+2:3" "2:+3" "1e0:2" "2:1e1" ":" "2:" ":3" " " "f:" ":L" "2:-0" "0x1:2" "f:f" "7:9" ///
    "8:8" "9:10" "1:7" "10" {
    local ++j
    pq_compare "IO-01 rowrange() form `j'" `"using "`fx'/rowgroups.parquet", clear rowrange(`r')"'
}
pq_compare "IO-01 rowrange() abbreviation rowr()" `"using "`fx'/rowgroups.parquet", clear rowr(2:9)"'
pq_compare "IO-01 rowrange() dictionary over row groups" `"using "`fx'/arrow_rowgroups.parquet", clear rowrange(4:6)"'
pq_compare "IO-01 rowrange() multi" `"using "`fx'/multi.parquet", clear rowrange(1050:2300)"'
pq_compare "IO-01 favormemory" `"using "`fx'/multi.parquet", clear favormemory"'
pq_compare "IO-01 favormemory with rowrange()" `"using "`fx'/rowgroups.parquet", clear favormemory rowrange(2:8)"'
pq_compare "IO-01 favormemory with columns" `"cat_unused plain_int cat using "`fx'/arrow.parquet", clear favormemory"'
pq_compare "IO-01 columnlist order" `"t s id using "`fx'/multi.parquet", clear"'
pq_compare "IO-01 columnlist repeated and unknown names" `"s zz id s using "`fx'/multi.parquet", clear"'
pq_compare "IO-01 columnlist none valid" `"zz yy using "`fx'/multi.parquet", clear"'
pq_compare "IO-01 columnlist unsupported only" `"st lst using "`fx'/nested.parquet", clear"'
pq_compare "IO-01 columnlist with unsupported" `"rep id m using "`fx'/nested.parquet", clear"'
pq_compare "IO-01 columnlist names v#" `"ok 1 u using "`fx'/names_select.parquet", clear"'
pq_compare "IO-01 columnlist duration" `"duration large_text using "`fx'/arrow.parquet", clear"'
pq_compare "IO-01 describe abbreviation" `"using "`fx'/types.parquet", desc"'
pq_compare "IO-01 describe with clear" `"using "`fx'/types.parquet", describe clear"'
pq_compare "IO-01 describe with columns" `"t_int8 using "`fx'/types.parquet", describe"'
pq_compare "IO-01 filename without using" `""`fx'/multi.parquet", clear"'
pq_compare "IO-01 filename without extension" `"using "`fx'/multi", clear"'
pq_compare "IO-01 describe without extension or using" `""`fx'/rowgroups", describe"'
pq_compare "IO-01 file not found" `"using "`fx'/no_such_file", clear"'
pq_compare "IO-01 file not found with extension" `"using "`fx'/no_such_file.parquet", clear"'
pq_compare "IO-01 URL" `"using "https://www.stata.com/examples/auto.parquet", clear"'
pq_compare "IO-01 option not allowed" `"using "`fx'/multi.parquet", clear bogus"'
capture program drop _pq_changed
program define _pq_changed
    sysuse auto, clear
    generate z = 1
end
pq_compare "IO-01 changed data need clear" `"using "`fx'/multi.parquet""' _pq_changed
pq_compare "IO-01 unchanged data need no clear" `"using "`fx'/multi.parquet""' `"sysuse auto, clear"'
pq_compare "IO-01 empty dataset needs no clear" `"using "`fx'/multi.parquet""' clear
* A path with spaces
local tmpdir `"`c(tmpdir)'/ctools parquet `=runiformint(1, 1e6)'"'
capture mkdir `"`tmpdir'"'
copy "`fx'/types.parquet" `"`tmpdir'/with space.parquet"', replace
pq_compare "IO-01 path with spaces" `"using "`tmpdir'/with space", clear"'
erase `"`tmpdir'/with space.parquet"'
capture rmdir `"`tmpdir'"'

* ---------------------------------------------------------------------------
* Brotli: pages from libbrotli at qualities 0-11 (brotli_types: every physical
* type, dictionary pages, data pages v1 and v2, row groups), files written by
* pyarrow, and stored meta-blocks (codec_brotli) are compared in the fixture
* loop above; here options over them, and damaged Brotli pages
* ---------------------------------------------------------------------------
foreach f in brotli_types brotli_pyarrow_v1 brotli_pyarrow_v2 corrupt_brotli {
    capture confirm file "`fx'/`f'.parquet"
    if _rc test_fail "IO-01 brotli fixture `f'" "missing: parquet_fixtures.py writes it when pyarrow is installed"
}
pq_compare "IO-01 brotli rowrange() across row groups" `"using "`fx'/brotli_types.parquet", clear rowrange(850:1950)"'
pq_compare "IO-01 brotli columnlist" `"b_text b_dict_text b_int96 b_int32_v2 using "`fx'/brotli_types.parquet", clear"'
pq_compare "IO-01 brotli favormemory" `"using "`fx'/brotli_types.parquet", clear favormemory"'
pq_compare "IO-01 brotli pyarrow rowrange()" `"using "`fx'/brotli_pyarrow_v2.parquet", clear rowrange(990:2010)"'
foreach how in "" snapshot {
    capture noisily {
        sysuse auto, clear
        tempfile before
        quietly save `before'
        capture noisily cimport parquet using "`fx'/corrupt_brotli.parquet", clear `how'
        assert _rc == 691
        cf _all using `before', all
        assert c(k) == 12 & _N == 74
    }
    if _rc test_fail "IO-01 damaged brotli fails cleanly, data kept `how'" "r(`=_rc')"
    else test_pass "IO-01 damaged brotli fails cleanly, data kept `how'"
}

* ---------------------------------------------------------------------------
* The undo copy (cimport's snapshot and fast options)
* ---------------------------------------------------------------------------
foreach how in "" snapshot {
    capture noisily {
        sysuse auto, clear
        tempfile before
        quietly save `before'
        capture noisily cimport parquet using "`fx'/names_reserved.parquet", clear `how'
        assert _rc == 3598
        cf _all using `before', all
        assert c(k) == 12
        capture noisily cimport parquet using "`fx'/empty.parquet", clear `how'
        assert _rc == 3900
        cf _all using `before', all
        global CTOOLS_cimport_fault load
        capture noisily cimport parquet using "`fx'/multi.parquet", clear `how'
        local rc = _rc
        global CTOOLS_cimport_fault
        assert `rc' == 1
        cf _all using `before', all
    }
    if _rc test_fail "IO-01 failure restores the data `how'" "r(`=_rc')"
    else test_pass "IO-01 failure restores the data `how'"
}
capture noisily {
    sysuse auto, clear
    capture noisily cimport parquet using "`fx'/names_reserved.parquet", clear fast
    assert _rc == 3598
    * fast keeps no copy: the failed import has replaced the data
    assert c(k) != 12
    sysuse auto, clear
    cimport parquet using "`fx'/multi.parquet", clear fast
    assert _N == 3000 & c(k) == 4
    capture noisily cimport parquet using "`fx'/multi.parquet", clear fast snapshot
    assert _rc == 198
}
if _rc test_fail "IO-01 fast keeps no copy; fast and snapshot conflict" "r(`=_rc')"
else test_pass "IO-01 fast keeps no copy; fast and snapshot conflict"
capture noisily {
    sysuse auto, clear
    cimport parquet using "`fx'/types.parquet", clear snapshot
    assert r(N) == 4 & r(k) == 42
}
if _rc test_fail "IO-01 snapshot on success" "r(`=_rc')"
else test_pass "IO-01 snapshot on success"

* ---------------------------------------------------------------------------
* IO-02/IO-03: FRED and Haver go to Stata's own commands, arguments unchanged.
* Only the paths that fail before any network or DLX access are exercised.
* ---------------------------------------------------------------------------
capture program drop _pq_delegated
program define _pq_delegated
    args label arguments
    tempfile reflog newlog
    sysuse auto, clear
    _pq_run `"import `arguments'"' `"`reflog'"'
    local rc1 = r(rc)
    local n1 = _N
    sysuse auto, clear
    _pq_run `"cimport `arguments'"' `"`newlog'"'
    local rc2 = r(rc)
    local n2 = _N
    mata: st_local("same", strofreal(_pq_same_output(st_local("reflog"), st_local("newlog"))))
    if `rc1' == `rc2' & `n1' == `n2' & `same' test_pass `"`label'"'
    else test_fail `"`label'"' "rc `rc1' vs `rc2', N `n1' vs `n2', same output `same'"
end
_pq_delegated "IO-02 fred no series" "fred"
_pq_delegated "IO-02 fred realtime() with vintage()" "fred GDP, realtime(2000-01-01 2001-01-01) vintage(2001-01-01) clear"
_pq_delegated "IO-02 fred realtime() with initial" "fred GDP, realtime(2000-01-01 2001-01-01) initial clear"
_pq_delegated "IO-02 fred allobs alone" "fred GDP, allobs clear"
_pq_delegated "IO-02 fred data in memory" "fred GDP UNRATE"
_pq_delegated "IO-02 fred serieslist() file" "fred, serieslist(ctools_no_such_file.txt) clear"
_pq_delegated "IO-02 fred daterange()" "fred GDP, daterange(notadate) clear"
_pq_delegated "IO-02 fred unknown option" "fred GDP, bogusoption clear"
_pq_delegated "IO-02 fred aggregate() with initial" "fred GDP, aggregate(annual) initial clear"
_pq_delegated "IO-03 haver no series" "haver"
_pq_delegated "IO-03 haver abbreviation" "hav"
_pq_delegated "IO-03 haver missing database" "haver gdp@ctools_no_such_db, clear"
_pq_delegated "IO-03 haver describe missing database" "haver gdp@ctools_no_such_db, describe"
_pq_delegated "IO-03 haver frommemory" "haver, frommemory"
_pq_delegated "IO-03 haver unknown option" "haver gdp@db, bogus"
_pq_delegated "IO-03 haverdirect no series" "haverdirect"
_pq_delegated "IO-03 haverdirect abbreviation" "haverd"
_pq_delegated "IO-03 haverdirect unknown option" "haverdirect gdp@usecon, bogus"

print_summary s2_parquet
if $TESTS_FAILED > 0 exit 9
global CTOOLS_COMPONENT_COMPLETE "s2_parquet"
    global CTOOLS_CASE_COMPLETE "s2_parquet"
    exit
}

di as error "Unknown s2_parquet case: `__case'"
exit 198
