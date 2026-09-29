do "validation/io_test_helpers.do" `0'

clear
set obs 5
gen byte small=_n
gen double precise=1.1234567890123456+_n
replace precise=.a in 4
replace precise=.z in 5
gen str12 text="abc"
replace text="" in 2
replace text=" abc " in 3
gen int day=td(01jan2020)+_n
format day %td
gen double stamp=clock("01jan2020 12:34:56","DMYhms")+_n
format stamp %tc
label define yesno 1 "Yes" 2 "No" 3 "Other"
label values small yesno
label variable precise "Precision label"
label data "Test data label"
save "temp/io_allformats/source.dta", replace

foreach pair in "spss sav" "sasxport8 v8xpt" {
    gettoken kind ext : pair
    local ext = strtrim("`ext'")
    export `kind' * using "temp/io_allformats/native/fixture.`ext'", replace
    cio_import_test using "temp/io_allformats/native/fixture.`ext'", kind(`kind') name("`kind' reader metadata and values")
    cio_import_test using "temp/io_allformats/native/fixture.`ext'", kind(`kind') name("`kind' upper case") opts(case(upper))
    use "temp/io_allformats/source.dta", clear
    cio_export_test, kind(`kind') ext(`ext') name("`kind' writer metadata and values")
}
cio_import_test using "temp/io_allformats/native/fixture.sav", kind(spss) name("SPSS selection before case and compression") opts(case(upper)) select(small precise if small>2 in 2/5)
use "temp/io_allformats/source.dta", clear
cio_export_test, kind(spss) ext(sav) name("SPSS no value labels") opts(novallabels)
cio_export_test, kind(spss) ext(sav) name("SPSS subset including missing") select(small precise text if small>2)
replace text="café Ω" in 1
label define yesno 1 "Oui café" 2 "Non Ω" 3 "Autre", replace
label variable text "Texte Ω"
label data "Données Ω"
foreach pair in "spss sav" "sasxport8 v8xpt" {
    gettoken kind ext : pair
    local ext=strtrim("`ext'")
    cio_export_test, kind(`kind') ext(`ext') name("`kind' Unicode writer metadata and values")
    cio_import_test using "temp/io_allformats/native/output.`ext'", kind(`kind') name("`kind' Unicode reader metadata and values")
    use "temp/io_allformats/source.dta", clear
    replace text="café Ω" in 1
    label define yesno 1 "Oui café" 2 "Non Ω" 3 "Autre", replace
    label variable text "Texte Ω"
    label data "Données Ω"
}
clear
set obs 3
gen strL text=12000*"a"
replace text="short" in 2
replace text="" in 3
cio_export_test, kind(excel) ext(xls) name("XLS long string continuation writer")
cio_import_test using "temp/io_allformats/result/output.xls", kind(excel) name("XLS long string continuation reader")
clear
set obs 2
gen strL text=6000*"éΩ"
replace text="日本語" in 2
cio_export_test, kind(excel) ext(xls) name("XLS Unicode continuation writer")
cio_import_test using "temp/io_allformats/result/output.xls", kind(excel) name("XLS Unicode continuation reader")
cio_export_test, kind(excel) ext(xlsx) name("XLSX long Unicode writer")
cio_import_test using "temp/io_allformats/result/output.xlsx", kind(excel) name("XLSX long Unicode reader")
use "temp/io_allformats/source.dta", clear
clear
set obs 3
gen strL text=2100*"a"
replace text="short" in 2
replace text="" in 3
cio_export_test, kind(spss) ext(sav) name("SPSS long string writer")
cio_export_test, kind(sasxport8) ext(v8xpt) name("XPORT8 long string writer")

use "temp/io_allformats/source.dta", clear
replace text=strtrim(text)
export sasxport5 * using "temp/io_allformats/native/fixture.xpt", replace
cio_import_test using "temp/io_allformats/native/fixture.xpt", kind(sasxport5) name("XPORT5 format label companion")
cio_import_test using "temp/io_allformats/native/fixture.xpt", kind(sasxport5) name("XPORT5 ignoring labels") opts(novallabels)
use "temp/io_allformats/source.dta", clear
replace text=strtrim(text)
cio_export_test, kind(sasxport5) ext(xpt) name("XPORT5 writer format label companion")
cio_export_test, kind(sasxport5) ext(xpt) name("XPORT5 no label companion") opts(vallabfile(none)) iopts(novallabels)
rename precise mylongvariablename
rename stamp mylongvariabletwo
cio_export_test, kind(sasxport5) ext(xpt) name("XPORT5 renamed collisions") opts(rename vallabfile(both))

/* An XPORT library may contain several members. */
export sasxport5 * using "temp/io_allformats/native/first.xpt", replace rename vallabfile(none)
replace small=small+10
export sasxport5 small text using "temp/io_allformats/native/second.xpt" in 1/3, replace rename vallabfile(none)
mata: fh=fopen("temp/io_allformats/native/first.xpt","r"); one=fread(fh,100000); fclose(fh)
mata: fh=fopen("temp/io_allformats/native/second.xpt","r"); two=fread(fh,100000); fclose(fh)
capture erase "temp/io_allformats/native/multiple.xpt"
mata: fh=fopen("temp/io_allformats/native/multiple.xpt","w"); fwrite(fh,one+substr(two,241,.)); fclose(fh)
cio_import_test using "temp/io_allformats/native/multiple.xpt", kind(sasxport5) name("XPORT5 first member") opts(novallabels)
cio_import_test using "temp/io_allformats/native/multiple.xpt", kind(sasxport5) name("XPORT5 member selection") opts(member(SECOND) novallabels)
cio_import_test using "temp/io_allformats/native/multiple.xpt", kind(sasxport5) name("XPORT5 lowercase member selection") opts(member(second) novallabels)
foreach member in "" "member(second)" {
    import sasxport5 "temp/io_allformats/native/multiple.xpt", describe `member'
    local n=r(N)
    local k=r(k)
    local size=r(size)
    local members `"`r(members)'"'
    local nmembers=r(n_members)
    capture noisily cimport sasxport5 "temp/io_allformats/native/multiple.xpt", describe `member'
    local rc=_rc
    if !`rc' {
        if r(N)!=`n' | r(k)!=`k' | r(size)!=`size' local rc=9
        if "`member'"=="" & (r(n_members)!=`nmembers' | `"`r(members)'"'!=`"`members'"') local rc=9
    }
    cio_record "XPORT5 describe results `member'" `rc'
}
foreach opt in "" "case(upper)" "case(lower)" {
    cio_companion_state_test using "temp/io_allformats/native/fixture.v8xpt", companion("temp/io_allformats/native/fixture.v8xpt") opts(`opt')
}
use "temp/io_allformats/source.dta", clear
replace text="café Ω" in 1
cexport sasxport8 * using "temp/io_allformats/result/companion.v8xpt", replace
cio_companion_state_test using "temp/io_allformats/native/fixture.v8xpt", companion("temp/io_allformats/result/companion.v8xpt")

use "temp/io_allformats/source.dta", clear
cexport sas * using "temp/io_allformats/result/fixture.sas7bdat", replace
cio_import_test using "temp/io_allformats/result/fixture.sas7bdat", kind(sas) name("SAS dataset reader against native")
cio_import_test using "temp/io_allformats/result/fixture.sas7bdat", kind(sas) name("SAS selection before case and compression") opts(case(lower)) select(small precise if small>2)

cio_import_test using "temp/io_allformats/fixture.zsav", kind(spss) name("ZSAV compression and long strings; valid display format") opts(zsav) longstring
cio_import_test using "temp/io_allformats/long.sas7bdat", kind(sas) name("SAS long strings; valid display format") longstring
cio_import_test using "temp/io_allformats/long.sas7bdat", kind(sas) name("SAS catalogue definitions; valid display format") longstring opts(bcat("temp/io_allformats/labels.sas7bcat"))
cio_import_test using "temp/io_allformats/formats.sav", kind(spss) name("SPSS display formats and user-defined missing values") nativebadfmt
cio_import_test using "temp/io_allformats/formats.sas7bdat", kind(sas) name("SAS standard and custom format names")

use "temp/io_allformats/source.dta", clear
drop stamp
export dbase * using "temp/io_allformats/native/fixture.dbf", replace
cio_import_test using "temp/io_allformats/native/fixture.dbf", kind(dbase) name("dBase reader exact values and fields")
cio_import_test using "temp/io_allformats/native/fixture.dbf", kind(dbase) name("dBase upper case") opts(case(upper))
use "temp/io_allformats/source.dta", clear
drop stamp
cio_export_test, kind(dbase) ext(dbf) name("dBase writer raw precision and dates")
foreach fmt in %10.2f %10.2e %10.0g {
    format precise `fmt'
    cio_export_test, kind(dbase) ext(dbf) name("dBase datafmt `fmt'") opts(datafmt)
}


use "temp/io_allformats/source.dta", clear
foreach header in default variables varlabels {
    local first ""
    if "`header'"!="default" local first firstrow(`header')
    cio_export_test, kind(excel) ext(xls) name("XLS `header' writer") opts(`first')
}
export excel using "temp/io_allformats/native/fixture.xls", firstrow(variables) sheet("Data sheet") replace
export excel small text using "temp/io_allformats/native/fixture.xls", sheet("Other", modify) cell(C3)
cio_workbook_test using "temp/io_allformats/native/fixture.xls", name("XLS workbook description preserves data")
cio_import_test using "temp/io_allformats/native/fixture.xls", kind(excel) name("XLS dates and literal worksheet") opts(firstrow sheet("Data sheet"))
cio_import_test using "temp/io_allformats/native/fixture.xls", kind(excel) name("XLS cell range") opts(cellrange(B2:D4) sheet("Data sheet"))
cio_import_test using "temp/io_allformats/native/fixture.xls", kind(excel) name("XLS explicit columns") opts(sheet("Data sheet")) select(identifier=A note=C)
cio_import_test using "temp/io_allformats/native/fixture.xls", kind(excel) name("XLS implicit columns") opts(sheet("Data sheet")) select(identifier number)
cio_import_test using "temp/io_allformats/native/fixture.xls", kind(excel) name("XLS all strings") opts(allstring firstrow sheet("Data sheet"))

use "temp/io_allformats/source.dta", clear
export excel using "temp/io_allformats/native/fixture.xlsx", firstrow(variables) sheet("Data sheet") replace
export excel small text using "temp/io_allformats/native/fixture.xlsx", sheet("Other", modify) cell(C3)
cio_workbook_test using "temp/io_allformats/native/fixture.xlsx", name("XLSX workbook description preserves data")
cio_import_test using "temp/io_allformats/native/fixture.xlsx", kind(excel) name("XLSX dates and literal worksheet") opts(firstrow sheet("Data sheet"))
cio_import_test using "temp/io_allformats/native/fixture.xlsx", kind(excel) name("XLSX cell range") opts(cellrange(B2:D4) sheet("Data sheet"))
cio_import_test using "temp/io_allformats/native/fixture.xlsx", kind(excel) name("XLSX explicit columns") opts(sheet("Data sheet")) select(identifier=A note=C)
cio_import_test using "temp/io_allformats/native/fixture.xlsx", kind(excel) name("XLSX all strings") opts(allstring firstrow sheet("Data sheet"))
clear
set obs 3
gen double exactfloat=_n+0.25
gen double empty=.
gen strL text=2100*"a"
export excel using "temp/io_allformats/native/long.xlsx", firstrow(variables) replace
cio_import_test using "temp/io_allformats/native/long.xlsx", kind(excel) name("XLSX long strings and types") opts(firstrow)

foreach shape in shape0 shape1 shape3 shape5 shape8 shape11 shape13 shape15 shape18 shape21 shape23 shape25 shape28 shape31 multipart nullrecord {
    cio_import_test using "temp/io_allformats/`shape'.shp", kind(shp) name("SHP `shape' reader binary headers")
    import shp using "temp/io_allformats/`shape'.shp", clear
    if "`shape'"=="shape0" continue
    if inlist("`shape'","shape28","shape31") {
        capture noisily export shp "temp/io_allformats/native/failure.shp", replace
        local native_rc=_rc
        capture noisily cexport shp "temp/io_allformats/result/failure.shp", replace
        local rc=_rc
        cio_record "SHP `shape' native unsupported measure export" `=cond(`rc'==`native_rc',0,9)'
    }
    else cio_export_test, kind(shp) ext(shp) name("SHP `shape' writer") opts(shx)
}

/* Compare failure state with native: Excel retains data; other readers clear. */
foreach pair in "spss sav" "shp shp" "dbase dbf" "excel xls" {
    gettoken kind ext : pair
    local ext = strtrim("`ext'")
    cio_import_state_test using "temp/io_allformats/corrupt.`ext'", kind(`kind') state(changed) opts(clear)
}
foreach pair in "sas sas7bdat" "spss sav" "sasxport5 xpt" "sasxport8 v8xpt" "dbase dbf" "xls xls" "xlsx xlsx" "shp shp" {
    gettoken format ext : pair
    local ext=strtrim("`ext'")
    local kind `format'
    if inlist("`format'","xls","xlsx") local kind excel
    foreach corruption in invalid truncated {
        foreach state in changed empty {
            cio_import_state_test using "temp/io_allformats/`format'_`corruption'.`ext'", kind(`kind') state(`state') opts(clear)
        }
    }
    if inlist("`format'","sas","spss","sasxport5","sasxport8") {
        cio_import_state_test using "temp/io_allformats/`format'_longinvalid.`ext'", kind(`kind') state(changed) opts(clear)
    }
}
cio_import_state_test using "temp/io_allformats/result/fixture.sas7bdat", kind(sas) state(changed) opts(clear encoding(bogus))
cio_import_state_test using "temp/io_allformats/native/fixture.sav", kind(spss) state(changed) opts(clear encoding(bogus))
foreach corruption in invalid longinvalid {
    cio_import_state_test using "temp/io_allformats/sasxport5_`corruption'.xpt", kind(sasxport5) state(saved)
}
/* XPORT5 uniquely permits replacement of saved, unchanged data. */
foreach pair in "delimited temp/io_allformats/memory.csv" "excel temp/io_allformats/native/fixture.xlsx" "excel temp/io_allformats/native/fixture.xls" "spss temp/io_allformats/native/fixture.sav" "sas temp/io_allformats/result/fixture.sas7bdat" "sasxport5 temp/io_allformats/native/fixture.xpt" "sasxport8 temp/io_allformats/native/fixture.v8xpt" "dbase temp/io_allformats/native/fixture.dbf" "shp temp/io_allformats/shape1.shp" {
    gettoken kind path : pair
    local path=strtrim("`path'")
    foreach state in saved changed zeroobs zerovars empty {
        cio_import_state_test using "`path'", kind(`kind') state(`state')
    }
}
di "IO_FORMATS_PASSED=$IO_FORMATS_PASSED FAILED=$IO_FORMATS_FAILED"
if $IO_FORMATS_FAILED exit 9
