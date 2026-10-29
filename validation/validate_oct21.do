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
