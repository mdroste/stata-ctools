/* Exact data and metadata comparisons against installed Stata. Run from root.
   No renaming, forced native-only options, or numerical tolerances. */
clear all
set more off
set linesize 255
adopath ++ "build"
args pluginpath
if `"`pluginpath'"' != "" adopath ++ `"`pluginpath'"'
capture mkdir "temp"
capture mkdir "temp/io_parity"
global IO_PARITY_PASSED 0
global IO_PARITY_FAILED 0

mata:
real scalar io_samefile(string scalar a, string scalar b) {
    real scalar x, y, same
    string matrix p, q
    x = fopen(a, "r"); y = fopen(b, "r"); same = 1
    do {
        p = fread(x, 65536); q = fread(y, 65536)
        if (!rows(p) | !rows(q)) {
            same = rows(p) == rows(q)
            break
        }
        if (p != q) {
            same = 0
            break
        }
    } while (strlen(p))
    fclose(x); fclose(y)
    return(same)
}
string matrix io_snapshot() {
    string matrix values
    real scalar j
    values = J(st_nobs(), st_nvar(), "")
    for (j=1; j<=st_nvar(); j++) {
        if (st_isnumvar(j)) values[,j] = strofreal(st_data(.,j), "%21x")
        else values[,j] = st_sdata(.,j)
    }
    return(values)
}

end

program define io_import_case
    syntax using/, NAME(string) [OPTS(string asis) EXCEL]
    local kind "delimited"
    if "`excel'" != "" local kind "excel"
    clear
    capture noisily import `kind'  using `"`using'"', `opts'
    local native_rc = _rc
    if `native_rc' {
        di as error "REFERENCE FAILED: `name': `native_rc'"
        global IO_PARITY_FAILED = $IO_PARITY_FAILED + 1
        exit
    }
    unab names : _all
    local n = _N
    local k = c(k)
    local i = 0
    foreach v of local names {
        local ++i
        local type`i' : type `v'
        local fmt`i' : format `v'
        local label`i' : variable label `v'
    }
    mata: io_reference = io_snapshot()
    clear
    capture noisily cimport `kind' using `"`using'"', `opts'
    local rc = _rc
    if !`rc' {
        if _N == `n' & c(k) == `k' {
            mata: st_local("same", strofreal(all(io_reference :== io_snapshot())))
            if !`same' local rc = 9
        }
        unab got : _all
        if "`got'" != "`names'" | _N != `n' | c(k) != `k' local rc = 9
        local i = 0
        foreach v of local got {
            local ++i
            local ty : type `v'
            local fm : format `v'
            local la : variable label `v'
            if "`ty'" != "`type`i''" | "`fm'" != "`fmt`i''" | `"`macval(la)'"' != `"`macval(label`i')'"' {
                di as error "METADATA: `v': `ty' / `type`i''; `fm' / `fmt`i''; `macval(la)' / `macval(label`i')'"
                local rc = 9
            }
        }
    }
    if `rc' {
        di as error "FAIL: `name' (rc=`rc')"
        global IO_PARITY_FAILED = $IO_PARITY_FAILED + 1
    }
    else {
        di as text "PASS: `name'"
        global IO_PARITY_PASSED = $IO_PARITY_PASSED + 1
    }
end

program define io_export_case
    syntax, NAME(string) [OPTS(string asis)]
    quietly export delimited using "temp/io_parity/reference.csv", replace `opts'
    capture noisily cexport delimited using "temp/io_parity/result.csv", replace `opts'
    local rc = _rc
    if !`rc' {
        mata: st_local("identical", strofreal(io_samefile("temp/io_parity/reference.csv", "temp/io_parity/result.csv")))
        if !`identical' local rc = 9
    }
    if `rc' {
        di as error "FAIL: `name' (rc=`rc')"
        global IO_PARITY_FAILED = $IO_PARITY_FAILED + 1
    }
    else {
        di as text "PASS: `name'"
        global IO_PARITY_PASSED = $IO_PARITY_PASSED + 1
    }
end

io_import_case using "temp/io_parity/headers.csv", name("header names and labels") opts(encoding(utf-8))
io_import_case using "temp/io_parity/noheaders.csv", name("automatic header detection") opts(encoding(utf-8))
io_import_case using "temp/io_parity/numbers.csv", name("numeric types and extended missings") opts(encoding(utf-8))
io_import_case using "temp/io_parity/quoted.csv", name("quoted numerics and exact string widths") opts(encoding(utf-8))
io_import_case using "temp/io_parity/outlier.csv", name("unsampled numeric outlier") opts(encoding(utf-8))
io_import_case using "temp/io_parity/exponents.csv", name("fractional scientific notation") opts(encoding(utf-8))

io_import_case using "temp/io_parity/numbers.csv", name("asfloat preserves integer storage") opts(encoding(utf-8) asfloat)
io_import_case using "temp/io_parity/numbers.csv", name("asdouble preserves integer storage") opts(encoding(utf-8) asdouble)
io_import_case using "temp/io_parity/quoted.csv", name("stringcols all") opts(encoding(utf-8) stringcols(_all))
io_import_case using "temp/io_parity/quoted.csv", name("stripquotes no") opts(encoding(utf-8) stripquotes(no))
io_import_case using "temp/io_parity/whitespace.csv", name("whitespace-only fields and padded missings") opts(encoding(utf-8))
io_import_case using "temp/io_parity/digits.csv", name("leading zeros and small exponents") opts(encoding(utf-8))
io_import_case using "temp/io_parity/multichunk.csv", name("parallel parse of mixed fields") opts(encoding(utf-8))

clear
set obs 18
gen double d = .
local i = 0
foreach n in 0.1 0.1234567890123456 1.234567890123456 123.4567890123456 0.0001 0.00001 1e16 1e20 1e-20 -0.1 1000000000 3.141592653589793 . .a .z 1e100 1e-100 9.999999999999999 {
    local ++i
    quietly replace d = `n' in `i'
}
gen float f = d
gen str1000 s = "hello world"
quietly replace s = "" in 1
quietly replace s = " abc " in 2
quietly replace s = "a,b" in 3
quietly replace s = char(34) + "hi" + char(34) in 4
quietly replace s = 200 * (char(34) + "ab,") in 5
io_export_case, name("default CSV bytes, precision, long quotes")
io_export_case, name("forced quoting") opts(quote)
io_export_case, name("tab delimiter") opts(delimiter(tab))
format d %20.5f
io_export_case, name("numeric display formats") opts(datafmt)
io_export_case, name("numeric display formats and quote") opts(datafmt quote)

clear
set obs 20000
set seed 723019
gen double d = (runiform() - .5) * 10^mod(_n, 40) / 1e20
gen float f = d
format d f %24.17e
quietly export delimited using "temp/io_parity/numbers_exact.csv", datafmt replace
io_export_case, name("randomized magnitudes and rounding")

/* Native quotes only fields containing the delimiter or a double quote. */
clear
set obs 8
gen str20 s = ""
quietly replace s = "a" + char(10) + "b" in 1
quietly replace s = "a" + char(13) + char(10) + "b" in 2
quietly replace s = "a" + char(9) + "b" in 3
quietly replace s = " lead" in 4
quietly replace s = "trail " in 5
quietly replace s = "a,b" in 6
quietly replace s = char(34) + "q" + char(34) in 7
quietly replace s = "  " in 8
io_export_case, name("quoting: delimiter and quote only")
io_export_case, name("quoting: tab delimiter") opts(delimiter(tab))

/* Date/time variables use their display format by default; formatted values
   are never quoted. Unlabeled values of labeled variables fall back to it. */
clear
set obs 6
gen int dday = mdy(1,4,2000) + _n
quietly replace dday = . in 5
quietly replace dday = .a in 6
format dday %td
gen double dtc = clock("2000-01-04 13:45:07.123", "YMDhms") + _n * 1000
format dtc %tc
gen int dtm = ym(2000, _n)
format dtm %tm
gen int dtq = yq(2000, mod(_n,4) + 1)
format dtq %tq
gen int dty = 1998 + _n
format dty %ty
gen int dleft = dday
format dleft %-td
gen int dcomma = dday
format dcomma %tdMonth_dd,_CCYY
io_export_case, name("date formats by default")
io_export_case, name("date formats under quote") opts(quote)
io_export_case, name("date formats with datafmt") opts(datafmt)
gen int dlab = dday
label define io_dlab .a "not asked"
label values dlab io_dlab
format dlab %tdCCYY-NN-DD
gen byte lab = _n
quietly replace lab = .b in 5
label define io_lab 1 "one" 3 "three, 3"
label values lab io_lab
io_export_case, name("labeled values with date and missing fallback")
io_export_case, name("labeled values with nolabel") opts(nolabel)

/* Long rows that the 500-row size sample never sees (every 60th row here). */
clear
set obs 30000
gen long id = _n
gen str1000 s = cond(inrange(_n, 20001, 22000) & mod(_n - 1, 60), 1000 * "x", "short")
io_export_case, name("rows longer than the size sample")

/* Excel's default export has no first row; explicit headings use names/labels. */
clear
set obs 3
gen byte id = _n
gen str8 text = "row " + string(_n)
label variable id "Identifier label"
label variable text "Text label"
foreach options in default variables varlabels {
    local opt ""
    if "`options'" != "default" local opt "firstrow(`options')"
    tempfile source native
    quietly save `source'
    quietly export excel using "temp/io_parity/native.xlsx", replace `opt'
    capture noisily cexport excel using "temp/io_parity/custom.xlsx", replace `opt'
    local rc = _rc
    if !`rc' {
        quietly import excel using "temp/io_parity/native.xlsx", clear
        quietly save `native'
        quietly import excel using "temp/io_parity/custom.xlsx", clear
        capture cf _all using `native', all
        local rc = _rc
    }
    if `rc' {
        di as error "FAIL: Excel `options' (rc=`rc')"
        global IO_PARITY_FAILED = $IO_PARITY_FAILED + 1
    }
    else {
        di as text "PASS: Excel `options'"
        global IO_PARITY_PASSED = $IO_PARITY_PASSED + 1
    }
    use `source', clear
}
quietly export excel using "temp/io_parity/sheetnames.xlsx", replace sheet("Data sheet")
io_import_case using "temp/io_parity/sheetnames.xlsx", name("Excel literal worksheet name") opts(sheet("Data sheet")) excel

/* The plugin must stay loaded after clear all: its OpenMP runtime registers
   fork handlers, and Stata used to crash on the next shell command. */
sysuse auto, clear
quietly cexport delimited using "temp/io_parity/clearall.csv", replace
local passed $IO_PARITY_PASSED
local failed $IO_PARITY_FAILED
clear all
shell echo IO_PARITY_SHELL_OK
global IO_PARITY_PASSED = `passed' + 1
global IO_PARITY_FAILED = `failed'
di as text "PASS: shell after clear all"

di "IO_PARITY_PASSED=$IO_PARITY_PASSED FAILED=$IO_PARITY_FAILED"
if $IO_PARITY_FAILED exit 9
