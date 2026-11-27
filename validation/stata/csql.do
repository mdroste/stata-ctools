* csql: SQL through DuckDB. Named cases share this file; each do-call has its
* own local scope.
gettoken __case 0 : 0
macro shift
if "`__case'" == "" local __case "csql"

if "`__case'" == "csql" {
do "validation/support/validate_setup.do"

* The DuckDB library: environment CSQL_DUCKDB, then the installed one. The
* suite is skipped where there is neither (DuckDB is not part of the
* package; csql, install downloads it).
local __lib : environment CSQL_DUCKDB
if `"`__lib'"' != "" global CSQL_DUCKDB `"`__lib'"'
capture quietly csql, version
if _rc {
    test_skip "csql" "no DuckDB library: type csql, install, or set CSQL_DUCKDB to its file"
    print_summary "csql"
    global CTOOLS_COMPONENT_COMPLETE "csql"
    global CTOOLS_CASE_COMPLETE "csql"
    exit
}
local __lib `"`r(library)'"'

tempfile __base
local __dir `"`__base'_csql"'
capture mkdir `"`__dir'"'

* _csql_sig: everything a dataset carries besides its values, in CSQL_SIG.
capture program drop _csql_sig
program define _csql_sig
    local sig "N=`=_N' k=`=c(k)'"
    foreach v of varlist * {
        local sig `"`sig' | `v' `: type `v'' `: format `v'' [`: value label `v''] {`: variable label `v''}"'
    }
    foreach v of varlist * {
        local l : value label `v'
        if "`l'" != "" {
            capture label list `l'
            if _rc == 0 local sig `"`sig' L[`l':`r(k)':`r(min)':`r(max)']"'
        }
    }
    global CSQL_SIG `"`sig'"'
end

* _csql_same <name> <file>: the data in memory equal those in file, values
* and description.
capture program drop _csql_same
program define _csql_same
    args name file
    _csql_sig
    local s1 `"$CSQL_SIG"'
    preserve
    quietly use `"`file'"', clear
    _csql_sig
    local s2 `"$CSQL_SIG"'
    restore
    capture cf _all using `"`file'"', all
    local rc = _rc
    if `rc' test_fail "`name'" "values differ (cf r(`rc'))"
    else if `"`s1'"' != `"`s2'"' test_fail "`name'" `"description differs: `s1' vs `s2'"'
    else test_pass "`name'"
end

* =============================================================================
* The data in memory: every storage type, missing values, labels, formats
* =============================================================================
clear
set obs 300
set seed 20261126
gen byte b = cond(_n == 1, -127, cond(_n == 2, 100, floor(runiform() * 50)))
gen int i = cond(_n == 1, -32767, cond(_n == 2, 32740, floor(runiform() * 2000)))
gen long l = cond(_n == 1, -2147483647, cond(_n == 2, 2147483620, floor(runiform() * 1e6)))
gen float f = runiform()
gen double d = rnormal() * 1e10
gen str4 s = "s" + string(_n) if mod(_n, 3)
gen strL bigs = "L" * (mod(_n, 5) * 700) + string(_n)
replace b = . in 5
replace i = .a in 6
replace l = .z in 7
replace f = .b in 8
replace d = . in 9
label define bl 0 "zero" 1 "one" 2 "two"
label values b bl
label variable d "a double"
format d %15.2f
format l %12.0gc
gen long d_day = mdy(1, 1, 1999) + _n
format d_day %td
gen double t_clock = clock("2001-02-03 04:05:06", "YMDhms") + 1000 * _n
format t_clock %tc
tempfile __mem __memref
quietly save `__mem'
* the result: every missing value is . (SQL has one NULL); the types are
* already the narrowest that hold the values, as csql stores them
foreach v of varlist b i l f d {
    quietly replace `v' = . if missing(`v')
}
quietly save `__memref'

use `__mem', clear
csql, query("SELECT * FROM stata") clear
_csql_same "csql [clear] every storage type, missing values and labels from memory" `__memref'
capture assert i[6] == . & l[7] == . & f[8] == .
if !_rc test_pass "csql extended missing values read as NULL and come back as ."
else test_fail "csql extended missing values read as NULL and come back as ." "values differ"

use `__mem', clear
csql, query("SELECT * FROM data") clear
_csql_same "csql data is the data in memory without a using file" `__memref'

use `__mem', clear
_csql_sig
local before `"$CSQL_SIG"'
csql, query("SELECT count(*) AS n, count(b) AS nb, sum(l) AS sl FROM stata") rows(5)
local n = r(N)
_csql_sig
quietly count if !missing(b)
local nb = r(N)
if `"`before'"' == `"$CSQL_SIG"' & `n' == 1 test_pass "csql [query] lists rows and leaves the data in memory alone"
else test_fail "csql [query] lists rows and leaves the data in memory alone" "data changed or r(N) = `n'"

* missing values and strings: NULL only for numeric missing
quietly use `__mem', clear
quietly count if missing(b)
local mb = r(N)
quietly count if missing(s)
local es = r(N)
quietly csql, query("SELECT count(*) - count(b) AS mb, count(*) - count(s) AS ms, sum(CASE WHEN s = '' THEN 1 ELSE 0 END) AS es FROM stata") clear
if mb[1] == `mb' & ms[1] == 0 & es[1] == `es' test_pass "csql missing values read as NULL, empty strings as ''"
else test_fail "csql missing values read as NULL, empty strings as ''" "mb `=mb[1]' (`mb'), ms `=ms[1]', es `=es[1]' (`es')"

* dates: %td and %tc are DATE and TIMESTAMP; nodatetime keeps numbers
use `__mem', clear
csql, query("SELECT d_day, typeof(d_day) AS td, t_clock, typeof(t_clock) AS tc, year(d_day) AS y, d_day + 1 AS d1, t_clock + INTERVAL 1 SECOND AS t1 FROM stata") clear
capture assert td == "DATE" & tc == "TIMESTAMP" & y == year(d_day) & d1 == d_day + 1 & t1 == t_clock + 1000
local ok = !_rc & "`: format d1'" == "%td" & "`: format t1'" == "%tc"
use `__mem', clear
csql, query("SELECT d_day, typeof(d_day) AS td FROM stata") clear nodatetime
capture assert td == "INTEGER" | td == "FLOAT" | td == "DOUBLE"
if `ok' & !_rc test_pass "csql [nodatetime] %td and %tc read as DATE and TIMESTAMP unless nodatetime"
else test_fail "csql [nodatetime] %td and %tc read as DATE and TIMESTAMP unless nodatetime" "types or values differ"

* labels follow names; nolabel drops them
use `__mem', clear
csql, query("SELECT b, d, d AS e FROM stata") clear
local ok = "`: value label b'" == "bl" & `"`: variable label d'"' == "a double" & "`: format d'" == "%15.2f" & `"`: variable label e'"' == ""
use `__mem', clear
csql, query("SELECT b, d FROM stata") clear nolabel
if `ok' & "`: value label b'" == "" & `"`: variable label d'"' == "" test_pass "csql [nolabel] labels follow names unless nolabel"
else test_fail "csql [nolabel] labels follow names unless nolabel" "labels not carried as documented"

* strl(): str# up to the threshold, strL past it
use `__mem', clear
csql, query("SELECT s, bigs, 'abcdefghij' AS t FROM stata") clear
local ok = "`: type s'" == "str4" & "`: type bigs'" == "strL" & "`: type t'" == "str10"
use `__mem', clear
csql, query("SELECT s, 'abcdefghij' AS t FROM stata") clear strl(5)
if `ok' & "`: type s'" == "str4" & "`: type t'" == "strL" test_pass "csql [strl] strings past strl() bytes are strL"
else test_fail "csql [strl] strings past strl() bytes are strL" "types `: type s' `: type t'"

* SQL names ignore case: a variable whose name differs from an earlier one
* only in case reads as name_2 (and keeps its labels under that name)
clear
set obs 3
gen x = _n
gen X = -_n
label variable X "upper"
csql, query("SELECT * FROM stata") clear
local names
foreach v of varlist * {
    local names "`names' `v'"
}
local names `names'
if "`names'" == "x X_2" & X_2[2] == -2 & `"`: variable label X_2'"' == "upper" ///
    test_pass "csql variables that differ only in case read as name_2"
else test_fail "csql variables that differ only in case read as name_2" "`names'"

* storage types like compress
csql, query("SELECT 100 AS a, 101 AS b, 32740 AS c, 32741 AS d, 2147483620 AS e, 2147483621 AS f, 1.5::FLOAT AS g, 1.5 AS h, true AS k, NULL AS m") clear
local types
foreach v of varlist * {
    local types "`types' `: type `v''"
}
local types `types'
if "`types'" == "byte int int long long double float double byte byte" test_pass "csql result types are the narrowest that hold the values"
else test_fail "csql result types are the narrowest that hold the values" "`types'"

* ENUM: integer codes with a value label
csql, query("CREATE TYPE mood AS ENUM ('sad', 'ok', 'happy'); SELECT 'happy'::mood AS m UNION ALL SELECT 'sad'::mood UNION ALL SELECT NULL") clear
capture assert m[1] == 3 & m[2] == 1 & m[3] == .
if !_rc & "`: value label m'" == "m" & "`: label m 3'" == "happy" test_pass "csql ENUM columns become labeled codes"
else test_fail "csql ENUM columns become labeled codes" "codes or label differ"

* nested and other types are converted to text
csql, query("SELECT [1, 2] AS lst, {'a': 1} AS st, INTERVAL 2 DAY AS iv") clear
capture assert lst == "[1, 2]" & st == "{'a': 1}" & iv == "2 days"
if !_rc test_pass "csql lists, structs and intervals become text"
else test_fail "csql lists, structs and intervals become text" "`=lst[1]' `=st[1]' `=iv[1]'"

* values Stata cannot hold
csql, query("SELECT 'nan'::DOUBLE AS a, 'inf'::DOUBLE AS b, 1e308 AS c, 9007199254740993::BIGINT AS d") clear
capture assert a == . & b == . & c == . & d == 9007199254740992
if !_rc test_pass "csql NaN, infinities and huge doubles become missing"
else test_fail "csql NaN, infinities and huge doubles become missing" "values differ"

* names
csql, query(`"SELECT 1 AS "my col", 2 AS "2x", 3 AS "_all", 4 AS "x", 5 AS "x", count(*)"') clear
local names
foreach v of varlist * {
    local names "`names' `v'"
}
local names `names'
if "`names'" == "my_col _2x __all x x_1 count_star__" & `"`: variable label my_col'"' == "my col" & "`r(renamed)'" != "" ///
    test_pass "csql invalid and repeated names become valid Stata names"
else test_fail "csql invalid and repeated names become valid Stata names" "`names'"

* =============================================================================
* Files: .dta in place, Parquet, CSV, JSON
* =============================================================================
sysuse auto, clear
quietly save `"`__dir'/auto118.dta"', replace
quietly saveold `"`__dir'/auto117.dta"', version(13) replace
quietly saveold `"`__dir'/auto115.dta"', version(12) replace
quietly saveold `"`__dir'/auto114.dta"', version(11) replace
* csql stores the narrowest types, as compress does; a left-aligned string
* format follows the new width
quietly compress
format make %-17s
tempfile __auto
quietly save `__auto'
foreach f in auto118 auto117 auto115 auto114 {
    csql, query("SELECT * FROM '`__dir'/`f'.dta'") clear
    _csql_same "csql reads `f' through FROM 'file.dta'" `__auto'
}
csql using `"`__dir'/auto118.dta"', query("SELECT * FROM data") clear
_csql_same "csql using file.dta reads it as data" `__auto'

csql, query(`"SELECT make, "foreign" FROM read_dta('`__dir'/auto118.dta', decode := true) WHERE make LIKE 'Volvo%'"') clear
if "`: type foreign'" == "str7" & foreign[1] == "Foreign" test_pass "csql read_dta decode := true gives label texts"
else test_fail "csql read_dta decode := true gives label texts" "`: type foreign' `=foreign[1]'"

* saving(): .dta, .parquet, .csv, .json
sysuse auto, clear
csql, query("SELECT * FROM stata") saving(`"`__dir'/a.dta"', replace)
use `__auto', clear
_csql_same "csql [saving] [replace] writes a .dta file" `"`__dir'/a.dta"'
sysuse auto, clear
csql, query("SELECT * FROM stata") saving(`"`__dir'/a.parquet"', replace)
csql using `"`__dir'/a.parquet"', query("SELECT * FROM data") clear
capture cf _all using `__auto'
if !_rc test_pass "csql saving() parquet round trip"
else test_fail "csql saving() parquet round trip" "values differ"
sysuse auto, clear
csql, query("SELECT make, price FROM stata") saving(`"`__dir'/a.csv"', replace)
csql, query("SELECT make, price FROM stata") saving(`"`__dir'/a.json"', replace)
csql using `"`__dir'/a.csv"', query("SELECT count(*) AS n, sum(price) AS p FROM data") clear
local ok = n[1] == 74 & p[1] == 456229
csql using `"`__dir'/a.json"', query("SELECT count(*) AS n, sum(price) AS p FROM data") clear
if `ok' & n[1] == 74 & p[1] == 456229 test_pass "csql saving() csv and json round trip"
else test_fail "csql saving() csv and json round trip" "counts differ"
capture csql, query("SELECT 1 AS x") saving(`"`__dir'/a.csv"')
if _rc == 602 test_pass "csql saving() without replace refuses an existing file"
else test_fail "csql saving() without replace refuses an existing file" "r(`=_rc')"
sysuse auto, clear
csql, query("SELECT make, price FROM stata WHERE price > 10000") saving(`"`__dir'/b.dta"', replace) clear
if _N == 10 & "`c(filename)'" != "" test_pass "csql saving() with clear loads the saved file"
else test_fail "csql saving() with clear loads the saved file" "N = `=_N'"

* collections of files
forvalues y = 2001/2003 {
    csql, query("SELECT `y' AS year, range AS id FROM range(10)") saving(`"`__dir'/part_`y'.parquet"', replace)
}
csql, query("SELECT year, count(*) AS n FROM '`__dir'/part_*.parquet' GROUP BY ALL ORDER BY year") clear
capture assert _N == 3 & n == 10 & year == 2000 + _n
if !_rc test_pass "csql reads a pattern of Parquet files as one table"
else test_fail "csql reads a pattern of Parquet files as one table" "N = `=_N'"

* strL and binary strL, in memory and in a file
clear
set obs 500
gen long id = _n
gen strL big = "x" * (2500 + mod(_n, 9)) + string(_n)
gen strL bin = "a" + char(0) + "b" if mod(_n, 2)
quietly save `"`__dir'/strl.dta"', replace
quietly compress
tempfile __strlref
quietly save `__strlref'
use `"`__dir'/strl.dta"', clear
csql, query("SELECT id, big, bin FROM stata") clear
_csql_same "csql strL and binary strL from memory" `__strlref'
csql, query("SELECT * FROM '`__dir'/strl.dta'") clear
_csql_same "csql strL and binary strL from a .dta file" `__strlref'

* order: a filter keeps the order of the rows it reads, at any size
clear
set obs 300000
gen long id = _n
gen double u = runiform()
gen str8 c = "c" + string(mod(_n, 17))
quietly save `"`__dir'/ord.dta"', replace
csql, query("SELECT * FROM stata WHERE u < 0.5") clear
capture assert id > id[_n - 1] if _n > 1
local ok = !_rc
csql, query("SELECT * FROM '`__dir'/ord.dta' WHERE u < 0.5") clear
capture assert id > id[_n - 1] if _n > 1
tempfile __ord
quietly save `__ord'
quietly use if u < 0.5 using `"`__dir'/ord.dta"', clear
capture cf _all using `__ord'
if `ok' & !_rc test_pass "csql filters keep row order, as use if"
else test_fail "csql filters keep row order, as use if" "rows out of order"

* =============================================================================
* Scripts, database files, engine options
* =============================================================================
tempfile __sqlf
tempname __fh
file open `__fh' using `__sqlf', write text replace
file write `__fh' "-- a comment; with a semicolon" _n
file write `__fh' "CREATE TABLE t AS SELECT range AS x FROM range(10);" _n
file write `__fh' "INSERT INTO t VALUES (100);" _n
file write `__fh' "SELECT sum(x) AS s, 'a;b' AS semi FROM t;" _n
file close `__fh'
csql, file(`__sqlf') clear
if s[1] == 145 & semi[1] == "a;b" test_pass "csql [file] runs an SQL script from a file"
else test_fail "csql [file] runs an SQL script from a file" "s = `=s[1]'"

capture erase `"`__dir'/w.duckdb"'
csql, database(`"`__dir'/w.duckdb"') query("CREATE TABLE firms AS SELECT range AS id, range % 3 AS ind FROM range(30)")
local changed = r(changed)
csql, database(`"`__dir'/w.duckdb"') query("SELECT ind, count(*) AS n FROM firms GROUP BY ALL ORDER BY ind") clear
if _N == 3 & n[1] == 10 & `changed' == 30 test_pass "csql [database] keeps tables between calls"
else test_fail "csql [database] keeps tables between calls" "N = `=_N'"
capture csql, database(`"`__dir'/w.duckdb"') readonly query("DROP TABLE firms")
local rc = _rc
csql, database(`"`__dir'/w.duckdb"') readonly query("SELECT count(*) AS n FROM firms") clear
if `rc' == 198 & n[1] == 30 test_pass "csql [readonly] opens the database for reading only"
else test_fail "csql [readonly] opens the database for reading only" "r(`rc')"

* memory() forces DuckDB to spill to tempdir(); threads() gives the same result
csql, query("SELECT count(*) AS n, sum(c) AS s FROM (SELECT r % 2000000 AS k, count(*) AS c FROM (SELECT range AS r FROM range(6000000)) GROUP BY k)") ///
    memory(64MB) tempdir(`"`__dir'/spill"') clear
local ok = n[1] == 2000000 & s[1] == 6000000
csql, query("SELECT count(*) AS n, sum(c) AS s FROM (SELECT r % 2000000 AS k, count(*) AS c FROM (SELECT range AS r FROM range(6000000)) GROUP BY k)") threads(1) clear
if `ok' & n[1] == 2000000 & s[1] == 6000000 test_pass "csql [memory] [tempdir] [threads] large aggregations under a memory limit"
else test_fail "csql [memory] [tempdir] [threads] large aggregations under a memory limit" "n = `=n[1]'"

* library(), version, verbose
csql, version library(`"`__lib'"')
local v `"`r(version)'"'
tempfile __log
quietly log using `__log', text name(csqlv) replace
csql, query("SELECT 1 AS x") library(`"`__lib'"') verbose
quietly log close csqlv
tempname __lh
local saw 0
file open `__lh' using `__log', read text
file read `__lh' line
while r(eof) == 0 {
    if strpos(`"`line'"', "csql timing:") local saw 1
    file read `__lh' line
}
file close `__lh'
if `"`v'"' != "" & `saw' test_pass "csql [version] [library] [verbose] report the library and timings"
else test_fail "csql [version] [library] [verbose] report the library and timings" "version `v', timing `saw'"

* install: the archive must match the pinned checksum; a library file is copied
tempfile __junk
file open `__fh' using `__junk', write text replace
file write `__fh' "not a DuckDB release" _n
file close `__fh'
copy `__junk' `"`__dir'/duckdb-shared-libs-windows-amd64.tar.gz"', replace
copy `__junk' `"`__dir'/duckdb-shared-libs-osx-arm64.tar.gz"', replace
copy `__junk' `"`__dir'/duckdb-shared-libs-osx-amd64.tar.gz"', replace
copy `__junk' `"`__dir'/duckdb-shared-libs-linux-amd64.tar.gz"', replace
capture csql, install from(`"`__dir'"') library(`"`__dir'/lib.bin"') replace
local rc1 = _rc
capture confirm file `"`__dir'/lib.bin"'
local absent = _rc != 0
csql, install from(`"`__lib'"') library(`"`__dir'/lib.bin"') replace
csql, version library(`"`__dir'/lib.bin"')
if `rc1' == 610 & `absent' & `"`r(version)'"' == `"`v'"' test_pass "csql [install] [from] [replace] checks the archive and copies a library"
else test_fail "csql [install] [from] [replace] checks the archive and copies a library" "r(`rc1'), absent `absent'"

* =============================================================================
* Errors leave the data alone
* =============================================================================
use `__mem', clear
_csql_sig
local before `"$CSQL_SIG"'
capture csql, query("SELECT * FROM nosuchtable") clear
local rc1 = _rc
capture csql, query("SELECT 1; SELECT * FROM nosuchtable; SELECT 2") clear
local rc2 = _rc
capture csql using `"`__dir'/nosuchfile.parquet"', query("SELECT * FROM data") clear
local rc3 = _rc
capture csql, query("CREATE TABLE t AS SELECT 1") clear
local rc4 = _rc
capture csql, clear
local rc5 = _rc
_csql_sig
capture cf _all using `__mem'
if `rc1' == 198 & `rc2' == 198 & `rc3' == 601 & `rc4' == 198 & `rc5' == 198 & !_rc & `"`before'"' == `"$CSQL_SIG"' ///
    test_pass "csql errors leave the data in memory unchanged"
else test_fail "csql errors leave the data in memory unchanged" "r(`rc1') r(`rc2') r(`rc3') r(`rc4') r(`rc5')"

clear
csql, query("SELECT * FROM range(0) t(x)") clear
if _N == 0 & c(k) == 1 test_pass "csql an empty result loads as 0 observations"
else test_fail "csql an empty result loads as 0 observations" "N = `=_N'"

print_summary "csql"
global CTOOLS_COMPONENT_COMPLETE "csql"
global CTOOLS_CASE_COMPLETE "csql"
exit
}

if "`__case'" == "csql_extra" {
do "validation/support/validate_setup.do"

* Extra checks against the same computation done with Stata commands
* (collapse, merge, joinby, sort, egen, summarize) and against the
* documented type, name, label and error contracts. Skipped, as the csql
* case is, where no DuckDB library is found.
local __lib : environment CSQL_DUCKDB
if `"`__lib'"' != "" global CSQL_DUCKDB `"`__lib'"'
capture quietly csql, version
if _rc {
    test_skip "csql_extra" "no DuckDB library: type csql, install, or set CSQL_DUCKDB to its file"
    print_summary "csql_extra"
    global CTOOLS_COMPONENT_COMPLETE "csql_extra"
    global CTOOLS_CASE_COMPLETE "csql_extra"
    exit
}

tempfile __base
local __dir `"`__base'_csqlx"'
capture mkdir `"`__dir'"'

* _cx_ok name rc [detail]: pass when rc is 0
capture program drop _cx_ok
program define _cx_ok
    args name rc detail
    if `rc' == 0 test_pass `"`name'"'
    else test_fail `"`name'"' `"`detail'"'
end

* _cx_cmp name file keys [tol]: the data in memory hold the values of file
* (the same variables, the same rows once both are sorted by keys; no keys:
* the same order); numbers agree within reldif() tol. Storage types may
* differ (csql stores the narrowest type).
capture program drop _cx_cmp
program define _cx_cmp
    args name file keys tol
    if "`tol'" == "" local tol 0
    local why
    preserve
    quietly {
        if "`keys'" != "" sort `keys', stable
        unab mine : _all
        local n1 = _N
        tempfile a
        save `a'
        use `"`file'"', clear
        unab theirs : _all
        local n2 = _N
        if "`keys'" != "" sort `keys', stable
    }
    local miss : list theirs - mine
    local extra : list mine - theirs
    if `n1' != `n2' local why "N `n1' vs `n2'"
    else if `"`miss'`extra'"' != "" local why "variables differ: missing `miss', extra `extra'"
    else {
        quietly {
            local j 0
            foreach v of local theirs {
                local ++j
                rename `v' __r`j'
            }
            merge 1:1 _n using `a', nogenerate
            local j 0
            foreach v of local theirs {
                local ++j
                capture confirm string variable `v'
                local s1 = !_rc
                capture confirm string variable __r`j'
                if `s1' != !_rc {
                    local why "`why' `v'(string vs numeric)"
                    continue
                }
                if `s1' count if `v' != __r`j'
                else count if !(`v' == __r`j' | reldif(`v', __r`j') <= `tol')
                if r(N) local why "`why' `v'(`r(N)' rows)"
            }
        }
    }
    restore
    if `"`why'"' == "" test_pass `"`name'"'
    else test_fail `"`name'"' `"`why'"'
end

* =============================================================================
* Numbers SQL makes that Stata can hold
* =============================================================================
* negative HUGEINT (what sum() of integers returns) and wide DECIMAL
csql, query("SELECT -1::HUGEINT AS a, -12345::HUGEINT AS b, -2000000000000::HUGEINT AS c, 7::HUGEINT AS d") clear
capture assert a == -1 & b == -12345 & c == -2e12 & d == 7
_cx_ok "csql negative HUGEINT values load exactly" `=_rc' "a=`=a[1]' b=`=b[1]' c=`=c[1]' (want -1 -12345 -2e12)"

csql, query("SELECT -1.5::DECIMAL(38,2) AS a, -7::DECIMAL(38,0) AS b, 2.25::DECIMAL(38,2) AS c") clear
capture assert a == -1.5 & b == -7 & c == 2.25
_cx_ok "csql negative DECIMAL(38,s) values load exactly" `=_rc' "a=`=a[1]' b=`=b[1]' (want -1.5 -7)"

csql, query("SELECT -1234.567::DECIMAL(18,3) AS a, 98765.4321::DECIMAL(9,4) AS b, -3::DECIMAL(4,0) AS c") clear
capture assert a == -1234.567 & b == 98765.4321 & c == -3 & "`: type c'" == "byte"
_cx_ok "csql DECIMAL up to 18 digits loads exactly, scale 0 as an integer type" `=_rc' "a=`=a[1]' b=`=b[1]' c=`=c[1]' `: type c'"

* sum() of byte/int/long variables (HUGEINT in DuckDB) against collapse
clear
set seed 4101
set obs 2000
gen int g = mod(_n, 20)
gen byte b = floor(runiform() * 120) - 110
gen int i = floor(runiform() * 4000) - 3000
gen long l = floor(runiform() * 2e6) - 1.5e6
tempfile __ints __intsref
quietly save `__ints'
collapse (sum) b i l, by(g)
quietly save `__intsref'
use `__ints', clear
csql, query("SELECT g, sum(b) AS b, sum(i) AS i, sum(l) AS l FROM stata GROUP BY g") clear
_cx_cmp "csql sum() of negative integer variables matches collapse (sum)" `__intsref' "g"

* integer type boundaries on the negative side, as compress
csql, query("SELECT -127 AS a, -128 AS b, -32767 AS c, -32768 AS d, -2147483647 AS e, -2147483648 AS f") clear
local types
foreach v of varlist * {
    local types "`types' `: type `v''"
}
local types `types'
capture assert a == -127 & b == -128 & c == -32767 & d == -32768 & e == -2147483647 & f == -2147483648
_cx_ok "csql negative integer limits get byte/int/long/double as compress" ///
    `=cond(_rc, _rc, "`types'" != "byte int int long long double")' "`types'"

* FLOAT beyond Stata's float range: missing (or double), as for doubles; not an error
capture csql, query("SELECT 3e38::FLOAT AS f, 1 AS one") clear
local rc = _rc
local ok = 0
if !`rc' local ok = missing(f[1]) | reldif(f[1], 3e38) < 1e-6
_cx_ok "csql FLOAT values past Stata's float range load (as missing) instead of failing" ///
    `=!`ok'' "r(`rc')"

* doubles down to c(mindouble) fit in Stata, past it they are missing
csql, query("SELECT -8.98e307 AS x, -9e307 AS y, 1e-310 AS z, 5e-324 AS w") clear
capture assert x == -8.98e307 & y == . & z == 1e-310 & w == 5e-324
_cx_ok "csql doubles near c(mindouble) and subnormals load; past it, missing" `=_rc' "x=`=x[1]' y=`=y[1]'"
clear
set obs 3
gen double x = -8.98e307 in 1
replace x = 5e-324 in 2
replace x = 8.98e307 in 3
tempfile __bigd
quietly save `__bigd'
csql, query("SELECT * FROM stata") clear
_cx_cmp "csql SELECT * keeps doubles down to c(mindouble)" `__bigd' ""

* NaN and infinities of FLOAT, booleans with NULL
csql, query("SELECT 'nan'::FLOAT AS a, '-inf'::FLOAT AS b, 1.5::FLOAT AS c, true AS t, false AS f, NULL::BOOLEAN AS n") clear
capture assert a == . & b == . & c == 1.5 & t == 1 & f == 0 & n == . & "`: type t'" == "byte"
_cx_ok "csql FLOAT NaN/inf become missing; BOOLEAN NULL is missing" `=_rc' "a=`=a[1]' b=`=b[1]'"

* =============================================================================
* Dates and times
* =============================================================================
csql, query("SELECT TIME '12:34:56.789' AS tm, TIMESTAMP '1950-06-15 12:00:00.123' AS ts, '2001-02-03 04:05:06.789'::TIMESTAMP_MS AS tms, '2001-02-03 04:05:06.789'::TIMESTAMP_NS AS tns, '2001-02-03 04:05:06'::TIMESTAMP_S AS tss, DATE '1900-01-01' AS d0, 'infinity'::DATE AS di, '-infinity'::TIMESTAMP AS tinf") clear
local c0 = clock("2001-02-03 04:05:06.789", "YMDhms")
capture assert tm == 45296789 & ts == clock("1950-06-15 12:00:00.123", "YMDhms") & tms == `c0' & tns == `c0' ///
    & tss == clock("2001-02-03 04:05:06", "YMDhms") & d0 == mdy(1, 1, 1900) & di == . & tinf == .
local ok = !_rc & "`: format tm'" == "%tcHH:MM:SS" & "`: format tns'" == "%tc" & "`: format d0'" == "%td"
_cx_ok "csql TIME, TIMESTAMP_S/MS/NS, pre-1960 dates and infinities" `=!`ok'' "tm=`=tm[1]' tms=`=tms[1]'"

capture csql, query("SELECT TIMESTAMPTZ '2020-01-01 00:00:00+02' AS t") clear
local rc = _rc
local ok = 0
if !`rc' local ok = t[1] == clock("2019-12-31 22:00:00", "YMDhms")
_cx_ok "csql TIMESTAMPTZ loads in UTC" `=!`ok'' "r(`rc')"

* %tm reads as a number and keeps its format; nodatetime keeps %tc numbers
clear
set obs 5
gen int m = tm(2000m1) + _n
format m %tm
gen double t = clock("2010-05-06 07:08:09", "YMDhms") + 1000 * _n
format t %tc
replace t = . in 3
tempfile __tm
quietly save `__tm'
csql, query("SELECT m, m + 1 AS m1, typeof(m) AS ty FROM stata") clear
capture assert m == tm(2000m1) + _n & m1 == m + 1 & ty == "SMALLINT"
local ok = !_rc & "`: format m'" == "%tm"
_cx_ok "csql %tm reads as a number and keeps its format" `=!`ok''
use `__tm', clear
keep t
tempfile __tmt
quietly save `__tmt'
csql, query("SELECT t, typeof(t) AS ty FROM stata") clear nodatetime
capture assert ty == "DOUBLE"
local ok = !_rc & "`: format t'" == "%tc"
drop ty
_cx_cmp "csql [nodatetime] %tc values come back unchanged as numbers" `__tmt' "" 0
_cx_ok "csql [nodatetime] a %tc variable reads as DOUBLE and keeps its format" `=!`ok''

* =============================================================================
* Strings
* =============================================================================
clear
set obs 6
gen str2045 s = ""
replace s = "é" * 1022 + "x" in 1
replace s = `"a,b "quoted" ; semi"' in 2
replace s = "  lead and trail  " in 3
replace s = "Ölwechsel ж 中文" in 4
replace s = "" in 5
replace s = "x" * 2045 in 6
gen long id = _n
tempfile __str
quietly save `__str'
csql, query("SELECT * FROM stata") clear
local ok = "`: type s'" == "str2045"
_cx_cmp "csql str2045, unicode, quotes and spaces round trip" `__str' ""
_cx_ok "csql a column of 2045-byte strings stays str2045" `=!`ok'' "`: type s'"

use `__str', clear
gen blen = strlen(s)
gen clen = ustrlen(s)
gen str2045 up = ustrupper(s)
tempfile __strref
quietly save `__strref'
use `__str', clear
csql, query("SELECT s, id, strlen(s) AS blen, length(s) AS clen, upper(s) AS up FROM stata") clear
_cx_cmp "csql strlen/length/upper on unicode text match strlen/ustrlen/ustrupper" `__strref' ""

csql, query("SELECT repeat('x', 2045) AS a, repeat('y', 2046) AS b, repeat('é', 5) AS c") clear strl(10)
local t1 "`: type a' `: type b' `: type c'"
csql, query("SELECT repeat('x', 2045) AS a, repeat('y', 2046) AS b, repeat('é', 5) AS c") clear
local t2 "`: type a' `: type b' `: type c'"
capture assert strlen(a) == 2045 & strlen(b) == 2046 & c == "ééééé"
_cx_ok "csql [strl] str# up to 2045 (or strl()) bytes, strL beyond" ///
    `=cond(_rc, _rc, "`t1'|`t2'" != "strL strL str10|str2045 strL str10")' "`t1'|`t2'"

csql, query("SELECT * FROM (VALUES (NULL::VARCHAR, '', 'a'), ('zz', '', NULL)) t(n, e, a)") clear
capture assert n[1] == "" & n[2] == "zz" & e == "" & a[2] == "" & "`: type e'" == "str1"
_cx_ok "csql NULL text is empty; an all-empty column is str1" `=_rc'

* =============================================================================
* Names
* =============================================================================
clear
set obs 4
gen byte order = _n
gen byte group = 10 * _n
gen byte select = -_n
gen byte from = 3
gen int abcdefghijklmnopqrstuvwxyz012345 = 100 + _n
gen int ABCDEFGHIJKLMNOPQRSTUVWXYZ012345 = 200 + _n
gen byte é = 7
tempfile __kw
quietly save `__kw'
csql, query("SELECT * FROM stata") clear
_cx_cmp "csql SQL keywords, 32-character case twins and unicode as variable names round trip" `__kw' ""
use `__kw', clear
csql, query(`"SELECT "order" + "group" AS t, "select" * "from" AS u FROM stata"') clear
capture assert t == 11 * _n & u == -3 * _n
_cx_ok "csql quoted keyword names compute as in Stata" `=_rc'

clear
set obs 3
gen x = _n
gen X = 10 * _n
gen X_2 = 100 * _n
csql, query("SELECT * FROM stata") clear
local names
foreach v of varlist * {
    local names "`names' `v'"
}
local names `names'
capture assert x == _n & X_2 == 10 * _n
local rc = _rc
capture confirm variable X_2_2
if !_rc & !`rc' capture assert X_2_2 == 100 * _n
_cx_ok "csql a case twin and a real name_2 variable keep their values" `=cond(`rc', `rc', _rc)' "`names'"

csql, query(`"SELECT 1 AS "abcdefghijklmnopqrstuvwxyz0123456789", 2 AS "abcdefghijklmnopqrstuvwxyz012345ZZ", 3 AS "_N", 4 AS "str5", 5 AS "strL", 6 AS "a", 7 AS "A", 8 AS "x y""') clear
local ren `"`r(renamed)'"'
local names
foreach v of varlist * {
    local names "`names' `v'"
}
local names `names'
local ok = "`names'" == "abcdefghijklmnopqrstuvwxyz012345 abcdefghijklmnopqrstuvwxyz0123_1 __N _str5 _strL a A x_y"
local ok = `ok' & `"`: variable label abcdefghijklmnopqrstuvwxyz0123_1'"' == "abcdefghijklmnopqrstuvwxyz012345ZZ"
local ok = `ok' & `"`: variable label __N'"' == "_N" & abcdefghijklmnopqrstuvwxyz0123_1[1] == 2 & A[1] == 7
local ok = `ok' & `: list sizeof ren' == 6
_cx_ok "csql long, reserved and case-different result names" `=!`ok'' "`names' | `ren'"

* =============================================================================
* Labels and formats
* =============================================================================
sysuse auto, clear
preserve
recast double price
gen byte n = 1
collapse (sum) n (mean) price, by(foreign)
tempfile __autoref
quietly save `__autoref'
restore
csql, query(`"SELECT "foreign", count(*) AS n, avg(price) AS price FROM stata GROUP BY "foreign""') clear
local ok = "`: value label foreign'" == "origin" & `"`: variable label foreign'"' == "Car origin" & `"`: variable label price'"' == "Price"
_cx_cmp "csql GROUP BY a labeled variable matches collapse" `__autoref' "foreign" 1e-12
_cx_ok "csql GROUP BY keeps the value and variable labels of the key" `=!`ok''

* an ENUM column and an inherited value label with the same name
clear
set obs 2
gen byte v = _n
label define m 1 "one" 2 "two"
label values v m
csql, query("SELECT 'b'::ENUM('a', 'b', 'c') AS m, v FROM stata") clear
local ok = "`: label (v) 1'" == "one" & "`: label (v) 2'" == "two" & "`: label (m) 2'" == "b"
_cx_ok "csql an ENUM label does not replace an inherited value label of the same name" ///
    `=!`ok'' "v shows `: label (v) 1'/`: label (v) 2', m shows `: label (m) 2'"

* two sources with different value labels of the same name
clear
set obs 2
gen byte k = _n
gen byte w = _n
label define bl 1 "file one" 2 "file two"
label values w bl
quietly save `"`__dir'/lab.dta"', replace
clear
set obs 2
gen byte k = _n
gen byte b = _n
label define bl 1 "mem one" 2 "mem two"
label values b bl
csql, query("SELECT s.k, s.b, f.w FROM stata s JOIN '`__dir'/lab.dta' f USING (k) ORDER BY k") clear
local ok = "`: label (b) 1'" == "mem one" & "`: label (w) 1'" == "file one"
_cx_ok "csql same-named value labels from two sources keep their own texts" ///
    `=!`ok'' "b: `: label (b) 1', w: `: label (w) 1'"

* unicode and 80-character variable labels, display formats
clear
set obs 3
gen double x = _n / 3
label variable x "Größe ж 中文 — unicode label"
gen long y = 1000 * _n
label variable y "12345678901234567890123456789012345678901234567890123456789012345678901234567890"
format x %9.3fc
format y %12.0gc
csql, query("SELECT x, y, x * 2 AS x2 FROM stata") clear
local ok = `"`: variable label x'"' == "Größe ж 中文 — unicode label" & strlen(`"`: variable label y'"') == 80
local ok = `ok' & "`: format x'" == "%9.3fc" & "`: format y'" == "%12.0gc" & `"`: variable label x2'"' == ""
_cx_ok "csql unicode and 80-character variable labels and formats carry over" `=!`ok''

* =============================================================================
* GROUP BY, joins, windows and sorts against Stata
* =============================================================================
clear
set seed 4102
set obs 100000
gen long id = _n
gen int g = floor(runiform() * 500)
replace g = . if mod(id, 997) == 0
gen double x = rnormal() * 100 if runiform() > 0.05
gen double u = runiform()
gen byte v = floor(runiform() * 20)
gen str8 s = cond(mod(g, 37) == 0, "", "s" + string(mod(g, 37)))
replace s = "é" + s if mod(g, 5) == 0
tempfile __x
quietly save `__x'

collapse (sum) sx = x (mean) mx = x (min) lo = x (max) hi = x (count) nx = x (sd) sdx = x (median) md = x, by(g)
tempfile __xref
quietly save `__xref'
use `__x', clear
csql, query("SELECT g, sum(x) AS sx, avg(x) AS mx, min(x) AS lo, max(x) AS hi, count(x) AS nx, stddev_samp(x) AS sdx, median(x) AS md FROM stata GROUP BY g") clear
_cx_cmp "csql GROUP BY with a missing key matches collapse (sum mean min max count sd median)" `__xref' "g" 1e-10

use `__x', clear
gen long n = 1
collapse (sum) n (mean) mu = u, by(s)
tempfile __sref
quietly save `__sref'
use `__x', clear
csql, query("SELECT s, count(*) AS n, avg(u) AS mu FROM stata GROUP BY s") clear
_cx_cmp "csql GROUP BY a string key with empty and unicode values matches collapse" `__sref' "s" 1e-12

* LEFT JOIN with a .dta file against merge m:1
clear
set obs 400
gen int g = _n - 1 + 50 * (_n > 200)
gen byte ind = mod(g, 7)
gen str8 nm = "f" + string(g)
quietly save `"`__dir'/firms.dta"', replace
use `__x', clear
quietly merge m:1 g using `"`__dir'/firms.dta"', keep(master match) nogenerate
tempfile __mref
quietly save `__mref'
use `__x', clear
csql, query("SELECT t.*, f.ind, f.nm FROM stata t LEFT JOIN '`__dir'/firms.dta' f USING (g)") clear
_cx_cmp "csql LEFT JOIN with a .dta file matches merge m:1 keep(master match)" `__mref' "id"

* the same join through Parquet and CSV copies of the file
csql using `"`__dir'/firms.dta"', query("SELECT * FROM data") saving(`"`__dir'/firms.parquet"', replace)
csql using `"`__dir'/firms.dta"', query("SELECT * FROM data") saving(`"`__dir'/firms.csv"', replace)
use `__x', clear
quietly merge m:1 g using `"`__dir'/firms.dta"', keep(match) nogenerate
keep id g ind nm
tempfile __m3
quietly save `__m3'
use `__x', clear
csql, query("SELECT t.id, t.g, p.ind, c.nm FROM stata t JOIN '`__dir'/firms.parquet' p USING (g) JOIN '`__dir'/firms.csv' c USING (g)") clear
_cx_cmp "csql a join over the data in memory, a Parquet and a CSV file matches merge" `__m3' "id"

* many-to-many INNER JOIN against joinby
clear
set seed 4103
set obs 300
gen long ida = _n
gen int k = floor(runiform() * 40)
gen double a = runiform()
tempfile __ja
quietly save `__ja'
clear
set obs 200
gen long idb = _n
gen int k = floor(runiform() * 45)
gen double bb = runiform()
quietly save `"`__dir'/jb.dta"', replace
use `__ja', clear
joinby k using `"`__dir'/jb.dta"'
tempfile __jref
quietly save `__jref'
use `__ja', clear
csql, query("SELECT t.*, u.idb, u.bb FROM stata t JOIN '`__dir'/jb.dta' u USING (k)") clear
_cx_cmp "csql many-to-many JOIN matches joinby" `__jref' "ida idb"

* FULL JOIN against merge 1:1
clear
set obs 300
gen long id = 2 * _n
gen double p = _n / 7
quietly save `"`__dir'/fp.dta"', replace
clear
set obs 250
gen long id = 3 * _n
gen str5 q = "q" + string(_n)
quietly save `"`__dir'/fq.dta"', replace
use `"`__dir'/fp.dta"', clear
quietly merge 1:1 id using `"`__dir'/fq.dta"', nogenerate
tempfile __fref
quietly save `__fref'
use `"`__dir'/fp.dta"', clear
csql, query("SELECT * FROM stata FULL JOIN '`__dir'/fq.dta' USING (id)") clear
_cx_cmp "csql FULL JOIN matches merge 1:1 keeping every row" `__fref' "id"

* window functions against by-group Stata
use `__x', clear
bysort g (id): gen long rn = _n
by g: gen double rs = sum(u)
by g: gen double lagu = u[_n - 1]
egen double tot = total(x), by(g) missing
egen long rk = rank(v), by(g) track
keep id g rn rs lagu tot rk
tempfile __wref
quietly save `__wref'
use `__x', clear
csql, query("SELECT id, g, row_number() OVER w AS rn, sum(u) OVER (w ROWS UNBOUNDED PRECEDING) AS rs, lag(u) OVER w AS lagu, sum(x) OVER (PARTITION BY g) AS tot, rank() OVER (PARTITION BY g ORDER BY v) AS rk FROM stata WINDOW w AS (PARTITION BY g ORDER BY id)") clear
_cx_cmp "csql row_number, running sum, lag, partition sum and rank match by-group Stata and egen" `__wref' "id" 1e-11

* ORDER BY against sort/gsort: binary string order, missing last
clear
set seed 4104
set obs 5000
gen long id = _n
gen str6 w = word("apple Apple banana Banana _x 9z é Zed zed", 1 + floor(runiform() * 9))
gen int v = floor(runiform() * 30) if runiform() > 0.1
tempfile __o
quietly save `__o'
gsort w -v id
tempfile __oref
quietly save `__oref'
use `__o', clear
csql, query("SELECT * FROM stata ORDER BY w, v DESC, id") clear
_cx_cmp "csql ORDER BY text, number DESC matches gsort (byte order, missing last)" `__oref' ""
use `__o', clear
sort v w id
quietly save `__oref', replace
use `__o', clear
csql, query("SELECT * FROM stata ORDER BY v, w, id") clear
_cx_cmp "csql ORDER BY number with NULLs, text matches sort" `__oref' ""

* count(DISTINCT) against egen tag
use `__x', clear
egen byte tg = tag(s g), missing
quietly count if tg
local nd = r(N)
csql, query("SELECT count(*) AS n FROM (SELECT DISTINCT s, g FROM stata)") clear
_cx_ok "csql SELECT DISTINCT over a missing key matches egen tag" `=n[1] != `nd'' "`=n[1]' vs `nd'"

* =============================================================================
* Row order, large N, threads
* =============================================================================
clear
set seed 4105
set obs 1000000
gen long id = _n
gen byte b = floor(runiform() * 100) if runiform() > 0.01
gen int i = floor(runiform() * 30000)
gen float f = runiform()
gen double d = rnormal() * 1e6
gen str10 s = "r" + string(floor(runiform() * 1e6))
gen long dt = mdy(1, 1, 1990) + floor(runiform() * 10000)
format dt %td
replace dt = . if mod(_n, 1000) == 0
tempfile __big
quietly save `__big'
csql, query("SELECT * FROM stata") clear
_cx_cmp "csql SELECT * of 1,000,000 rows keeps every value and the row order" `__big' ""

use `__big', clear
quietly keep if f < 0.3 & !missing(b)
tempfile __bigf
quietly save `__bigf'
use `__big', clear
csql, query("SELECT * FROM stata WHERE f < 0.3 AND b IS NOT NULL") clear
_cx_cmp "csql WHERE on 1,000,000 rows keeps the order of keep if" `__bigf' ""

use `__big', clear
quietly keep in 500001/500010
quietly save `__bigf', replace
use `__big', clear
csql, query("SELECT * FROM stata LIMIT 10 OFFSET 500000") clear
_cx_cmp "csql LIMIT/OFFSET without ORDER BY matches keep in" `__bigf' ""

use `__big', clear
gen long n = 1
collapse (sum) n (mean) md = d (max) mxs = i, by(b)
quietly save `__bigf', replace
use `__big', clear
csql, query("SELECT b, count(*) AS n, avg(d) AS md, max(i) AS mxs FROM stata GROUP BY b") clear
_cx_cmp "csql GROUP BY over 1,000,000 rows matches collapse" `__bigf' "b" 1e-11

use `__big', clear
csql, query("SELECT b, count(*) AS n, min(d) AS lo, max(s) AS ms FROM stata WHERE i > 100 GROUP BY b ORDER BY b") clear threads(1)
tempfile __t1
quietly save `__t1'
use `__big', clear
csql, query("SELECT b, count(*) AS n, min(d) AS lo, max(s) AS ms FROM stata WHERE i > 100 GROUP BY b ORDER BY b") clear threads(8)
capture cf _all using `__t1', all
_cx_ok "csql [threads] threads(1) and threads(8) give the same rows" `=_rc'

* =============================================================================
* Other sources: frames, files with missing values, CSV text, empty data
* =============================================================================
clear
set obs 2
gen a = 1
tempfile __fr0
quietly save `__fr0'
capture frame drop cxframe
frame create cxframe
frame cxframe {
    quietly set obs 3
    gen z = 10 * _n
}
frame cxframe: csql, query("SELECT sum(z) AS s, count(*) AS n FROM stata") clear
frame cxframe: scalar __cx_s = s[1]
local s = scalar(__cx_s)
scalar drop __cx_s
capture cf _all using `__fr0', all
_cx_ok "csql reads the current frame and leaves the others alone" `=cond(_rc, _rc, `s' != 60)' "s = `s'"
frame drop cxframe

clear
set obs 30
gen double y = _n
replace y = .a in 1/3
replace y = .z in 4
replace y = . in 5
gen str3 t = cond(_n <= 2, "", "t")
quietly save `"`__dir'/mv.dta"', replace
csql, query("SELECT count(*) - count(y) AS ny, sum(CASE WHEN t = '' THEN 1 ELSE 0 END) AS ne FROM '`__dir'/mv.dta'") clear
_cx_ok "csql a .dta file's .a-.z read as NULL and empty strings as ''" `=cond(ny[1] != 5, 1, ne[1] != 2)' "ny=`=ny[1]' ne=`=ne[1]'"

clear
set obs 5
gen str20 t = ""
replace t = "a,b" in 1
replace t = `"say "hi""' in 2
replace t = "  lead" in 3
replace t = "trail  " in 4
replace t = "é ü ж" in 5
gen double num = _n / 4
tempfile __csv
quietly save `__csv'
csql, query("SELECT * FROM stata") saving(`"`__dir'/t.csv"', replace)
csql using `"`__dir'/t.csv"', query("SELECT * FROM data") clear
_cx_cmp "csql saving() .csv and using .csv round trip commas, quotes and spaces" `__csv' ""

clear
set obs 4
gen long id = _n
gen str5 nm = "n" + string(_n)
quietly save `"`__dir'/u.dta"', replace
clear
set obs 3
gen z = 1
csql using `"`__dir'/u.dta"', query("SELECT (SELECT count(*) FROM data) AS nd, (SELECT count(*) FROM stata) AS ns") clear
_cx_ok "csql using: data is the file and stata the data in memory in one query" `=cond(nd[1] != 4, 1, ns[1] != 3)'

* 0 observations: in memory, in a file, in a result
use `"`__dir'/u.dta"', clear
quietly drop in 1/l
quietly save `"`__dir'/empty.dta"', replace
csql, query("SELECT * FROM stata") clear
local ok = _N == 0 & c(k) == 2
capture confirm string variable nm
local ok = `ok' & !_rc
csql, query("SELECT * FROM '`__dir'/empty.dta'") clear
local ok = `ok' & _N == 0 & c(k) == 2
capture confirm string variable nm
local ok = `ok' & !_rc
csql, query("SELECT 'abc' AS s, 1 AS x WHERE false") clear
local ok = `ok' & _N == 0 & c(k) == 2
capture confirm string variable s
_cx_ok "csql 0-observation data, files and results load with their columns" `=cond(_rc, _rc, !`ok')'

clear
quietly set obs 5
capture csql, query("SELECT * FROM stata") clear
_cx_ok "csql table stata with observations but no variables is an error that keeps _N" ///
    `=cond(_rc == 0, 1, _N != 5 | c(k) != 0)' "r(`=_rc') N=`=_N'"

clear
set obs 1
gen x = 42
csql, query("SELECT x, x + 1 AS y FROM stata") clear
_cx_ok "csql one observation" `=cond(_N != 1, 1, y[1] != 43)'

* listing: r(N) counts every row; rows(0) lists none
quietly csql, query("SELECT range AS r, 'x' AS s FROM range(1000)") rows(3)
local n1 = r(N)
local k1 = r(k)
quietly csql, query("SELECT range AS r FROM range(77)") rows(0)
_cx_ok "csql [rows] r(N) and r(k) count all rows of a listing" `=cond(`n1' != 1000 | `k1' != 2, 1, r(N) != 77)' "`n1' `k1' `r(N)'"

* =============================================================================
* Errors and saving()
* =============================================================================
use `__x', clear
quietly datasignature
local sig0 `"`r(datasignature)'"'
capture csql, query("SELECT error('boom') AS e") clear
local rc1 = _rc
capture csql, query("SELECT id, CAST(s || 'x' AS INTEGER) AS bad FROM stata") clear
local rc2 = _rc
capture csql, query("SELECT CASE WHEN id = 90000 THEN error('late') ELSE id END AS id FROM stata") clear
local rc3 = _rc
quietly datasignature
local sig1 `"`r(datasignature)'"'
capture cf _all using `__x', all
_cx_ok "csql runtime SQL errors (error(), casts, late in the scan) leave the data unchanged" ///
    `=cond(_rc, _rc, `rc1' != 198 | `rc2' != 198 | `rc3' != 198 | `"`sig1'"' != `"`sig0'"')' "r(`rc1') r(`rc2') r(`rc3')"

* a failing query does not damage an existing saving() file
csql, query("SELECT range AS r FROM range(50)") saving(`"`__dir'/keep.dta"', replace)
csql, query("SELECT range AS r FROM range(50)") saving(`"`__dir'/keep.parquet"', replace)
quietly checksum `"`__dir'/keep.dta"'
local ck1 = r(checksum)
quietly checksum `"`__dir'/keep.parquet"'
local ck2 = r(checksum)
capture csql, query("SELECT CASE WHEN id = 90000 THEN error('late') ELSE id END AS id FROM stata") saving(`"`__dir'/keep.dta"', replace)
local rc1 = _rc
capture csql, query("SELECT CASE WHEN id = 90000 THEN error('late') ELSE id END AS id FROM stata") saving(`"`__dir'/keep.parquet"', replace)
local rc2 = _rc
quietly checksum `"`__dir'/keep.dta"'
local same = r(checksum) == `ck1'
capture quietly checksum `"`__dir'/keep.parquet"'
local same = `same' & !_rc & r(checksum) == `ck2'
_cx_ok "csql [saving] a failing query leaves an existing .dta / .parquet file intact" ///
    `=cond(`rc1' == 0 | `rc2' == 0, 1, !`same')' "r(`rc1') r(`rc2') same `same'"

* option errors
capture csql, query("SELECT 1 AS x") strl(0)
local r1 = _rc
capture csql, query("SELECT 1 AS x") rows(-1)
local r2 = _rc
capture csql, query("SELECT 1 AS x") threads(-1)
local r3 = _rc
capture csql, query("SELECT 1 AS x") saving(`"`__dir'/bad.xls"')
local r4 = _rc
capture confirm file `"`__dir'/bad.xls"'
local r4b = _rc
capture csql, query("SELECT 1 AS x") saving(`"`__dir'/c.parquet"', replace) clear
local r5 = _rc
capture csql, query("SELECT 1 AS x") file(`"`__dir'/nosuch.sql"')
local r6 = _rc
capture csql, file(`"`__dir'/nosuch.sql"')
local r7 = _rc
capture csql, query("SELECT 1 AS x") readonly
local r8 = _rc
use `__x', clear
_cx_ok "csql bad options fail with the documented codes and keep the data" ///
    `=cond("`r1' `r2' `r3' `r4' `r5' `r6' `r7' `r8'" != "198 198 198 198 198 198 601 198" | !`r4b', 1, 0)' ///
    "`r1' `r2' `r3' `r4' (`r4b') `r5' `r6' `r7' `r8'"

* saving(): a name without extension is a .dta file; r() results
capture erase `"`__dir'/noext.dta"'
csql, query("SELECT range AS r, 'a' AS s FROM range(12)") saving(`"`__dir'/noext"')
local ok = r(N) == 12 & r(k) == 2 & `"`r(saving)'"' == `"`__dir'/noext.dta"'
capture confirm file `"`__dir'/noext.dta"'
_cx_ok "csql [saving] a name without extension saves .dta; r(N) r(k) r(saving)" `=cond(_rc, _rc, !`ok')' `"`r(saving)'"'

* parquet round trip keeps dates, strings and missing values
use in 1/5000 using `__big', clear
csql, query("SELECT * FROM stata") saving(`"`__dir'/big.parquet"', replace)
csql using `"`__dir'/big.parquet"', query("SELECT * FROM data") clear
local ok = "`: format dt'" == "%td"
tempfile __b5
preserve
use in 1/5000 using `__big', clear
quietly save `__b5'
restore
_cx_cmp "csql saving() .parquet round trip of every type with missing values" `__b5' ""
_cx_ok "csql a DATE read back from Parquet is %td" `=!`ok''

* a dollar sign written \$ reaches the SQL
csql, query("SELECT 'a\$b' AS s, length('a\$b') AS n") clear
_cx_ok "csql \$ in query() is a literal dollar sign" `=cond(s[1] != "a" + char(36) + "b", 1, n[1] != 3)' "`=s[1]'"

* file(): the SQL is not macro-expanded
tempfile __sqlx
tempname __fh
file open `__fh' using `__sqlx', write text replace
file write `__fh' "SELECT 'a" _char(36) "b" _char(96) "c' || chr(39) AS s;" _n
file close `__fh'
csql, file(`__sqlx') clear
_cx_ok "csql [file] SQL from a file keeps dollar signs and backquotes" ///
    `=s[1] != "a" + char(36) + "b" + char(96) + "c" + char(39)' "`=s[1]'"

* extended missing values are one NULL group
clear
set obs 40
gen byte k = mod(_n, 4)
replace k = .a in 1/3
replace k = .b in 4/6
replace k = . in 7
csql, query("SELECT count(*) AS n FROM stata WHERE k IS NULL") clear
_cx_ok "csql . and .a-.z are all NULL in WHERE ... IS NULL" `=n[1] != 7' "n = `=n[1]'"

* a table kept in database() round trips the data in memory
capture erase `"`__dir'/rt.duckdb"'
use in 1/20000 using `__big', clear
csql, database(`"`__dir'/rt.duckdb"') query("CREATE TABLE t AS SELECT * FROM stata")
local ch = r(changed)
tempfile __rt
quietly save `__rt'
clear
csql, database(`"`__dir'/rt.duckdb"') query("SELECT * FROM t") clear
_cx_cmp "csql [database] a table made from the data in memory reads back unchanged" `__rt' ""
_cx_ok "csql [database] r(changed) counts the rows CREATE TABLE AS stores" `=`ch' != 20000' "`ch'"

print_summary "csql_extra"
global CTOOLS_COMPONENT_COMPLETE "csql_extra"
global CTOOLS_CASE_COMPLETE "csql_extra"
exit
}

di as error "Unknown csql case: `__case'"
exit 198
