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

di as error "Unknown csql case: `__case'"
exit 198
