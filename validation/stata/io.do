* io: named cases share this file, each do-call has its own local scope.
gettoken __case 0 : 0
macro shift
if "`__case'" == "" local __case "io"

if "`__case'" == "io_delimited_options" {
do "validation/support/io_test_helpers.do"
args root
if `"`root'"' == "" local root "`c(pwd)'/temp/io_delimited_options"
foreach opt in `"delimiters("|;",collapse)"' `"delimiters("||",asstring)"' `"delimiters("|;")"' `"delimiters("|;") collapsedelimiters"' {
    cio_import_test using "`root'/multidelim.txt", kind(delimited) name("delimiter case") opts(varnames(nonames) `macval(opt)')
}
cio_import_test using "`root'/long.csv", kind(delimited) name("CSV long Unicode strL") opts(encoding(utf8))
cio_import_test using "`root'/long.csv", kind(delimited) name("CSV long Unicode favorstrfixed") opts(encoding(utf8) favorstrfixed)
foreach names in "one" "one two" {
    cio_import_test using "`root'/multiline.csv", kind(delimited) name("CSV explicit names `names'") select(`names')
}
foreach limit in 1 2 3 20 unlimited {
    capture noisily import delimited using "`root'/multiline.csv", bindquotes(strict) maxquotedrows(`limit') clear
    local nr=_rc
    capture noisily cimport delimited using "`root'/multiline.csv", bindquotes(strict) maxquotedrows(`limit') clear
    local cr=_rc
    local rc=cond(`nr'==`cr',0,9)
    cio_record "maxquotedrows(`limit') error" `rc'
}

foreach opt in "varnames(1) rowrange(3:4)" "varnames(1) rowrange(2:3)" "rowrange(3:4)" "varnames(1) rowrange(4)" "varnames(1) rowrange(1:1)" {
    cio_import_test using "`root'/range.csv", kind(delimited) name("CSV range inference `opt'") opts(`opt')
}
import delimited using "`root'/long.csv", encoding(utf8) clear
cio_export_test, kind(delimited) ext(csv) name("CSV long Unicode writer") opts(quote) iopts(encoding(utf8))
import delimited using "`root'/long.csv", encoding(utf8) clear
cio_export_test, kind(delimited) ext(csv) name("CSV long Unicode filtered writer") select(text if id==1) iopts(encoding(utf8))
foreach enc in utf-32le UTF-32BE utf-16be utf16 windows-1251 cp1251 SHIFT-JIS shift_jis sjis ISO-8859-2 latin2 {
    local f cp1251
    if strpos(lower("`enc'"),"32le") local f utf32le
    if strpos(lower("`enc'"),"32be") local f utf32be
    if strpos(lower("`enc'"),"16") local f utf16be
    if strpos(lower("`enc'"),"jis") | "`enc'"=="sjis" local f sjis
    if strpos(lower("`enc'"),"8859") | "`enc'"=="latin2" local f latin2
    cio_import_test using "`root'/`f'.csv", kind(delimited) name("CSV encoding `enc'") opts(encoding(`enc'))
}
cio_import_test using "`root'/utf32bom.csv", kind(delimited) name("CSV UTF-32 BOM detection")
foreach opt in "encoding(bogus)" "encoding(utf32le)" "encoding(utf16be)" "locale(bogus)" "charset(cp1252) encoding(utf8)" {
    capture noisily import delimited using "`root'/range.csv", `opt' clear
    local nr=_rc
    capture noisily cimport delimited using "`root'/range.csv", `opt' clear
    local cr=_rc
    cio_record "CSV option error `opt'" `=cond(`nr'==`cr',0,9)'
}
foreach width in 20 40 50 100 200 1000 {
    foreach pct in 1 28 29 34 35 50 63 64 65 66 67 70 90 93 94 95 {
        cio_import_test using "`root'/storagew`width'p`pct'.csv", kind(delimited) name("CSV storage width `width' fraction `pct'")
    }
    cio_import_test using "`root'/storagew`width'p1.csv", kind(delimited) name("CSV favorstrfixed width `width'") opts(favorstrfixed)
}
foreach pair in "us en_US" "german de" "german de_DE" "french fr_FR" "french_nbsp fr_FR" "russian ru_RU" "swiss de_CH" "arabic ar_EG" "arabic_minus ar_EG" "persian fa_IR" "persian_minus fa_IR" "hindi hi_IN" "bad_sign en_US" "scientific en_US" "infinity en_US" "bad_groups en_US" "groups en_US" "trailing_group en_US" {
    gettoken file locale : pair
    local locale=strtrim("`locale'")
    cio_import_test using "`root'/`file'.csv", kind(delimited) name("CSV locale `file' `locale'") opts(delimiters(";") parselocale(`locale') encoding(utf8))
}
foreach opt in "stripquotes(yes)" "stripquotes(no)" "stripquotes(default)" "bindquotes(nobind)" "rowrange(f:l)" "colrange(f:l)" {
    cio_import_test using "`root'/quotes.csv", kind(delimited) name("CSV quotes/ranges `opt'") opts(`opt')
}
foreach opt in "maxquotedrows(1.5)" "maxquotedrows(bogus)" "rowrange(0)" "rowrange(4:2)" "rowrange(:)" "colrange(0)" "colrange(4:2)" "parselocale(bogus)" "parselocale(en_XX)" {
    capture noisily import delimited using "`root'/range.csv", `opt' clear
    local nr=_rc
    capture noisily cimport delimited using "`root'/range.csv", `opt' clear
    local cr=_rc
    cio_record "CSV invalid range/locale `opt'" `=cond(`nr'==`cr',0,9)'
}
foreach opt in "stripquotes(yes)" "stripquotes(no)" "stripquotes(default)" "bindquotes(nobind)" "bindquotes(nobind) stripquotes(yes)" "numericcols(2) stripquotes(no)" {
    cio_import_test using "`root'/quote_numbers.csv", kind(delimited) name("CSV numeric quotes `opt'") opts(`opt')
}
foreach enc in utf8 UTF-8 windows-1252 cp1252 latin1 Latin1 LATIN1 latin2 latin9 ascii US-ASCII {
    cio_import_test using "`root'/range.csv", kind(delimited) name("CSV returned encoding `enc'") opts(encoding(`enc'))
}
foreach opt in "" "encoding(utf8)" "encoding(latin1)" {
    cio_import_test using "`root'/empty.csv", kind(delimited) name("CSV empty returned results `opt'") opts(`opt')
}
foreach opt in `"delimiters("tab")"' "delimiters(tab)" `"delimiters("comma")"' "delimiters(comma)" `"delimiters("space")"' "delimiters(space)" "delimiters(whitespace)" {
    cio_import_test using "`root'/words.txt", kind(delimited) name("CSV literal and shortcut delimiter") opts(`opt')
}
foreach kind in tab semi pipe colon space {
    cio_import_test using "`root'/autodelim_`kind'.csv", kind(delimited) name("CSV automatic delimiter `kind'")
}
foreach kind in western central cyrillic koi8 arabic greek hebrew turkish japanese eucjp korean chinese traditional iso2022 {
    cio_import_test using "`root'/detect_`kind'.csv", kind(delimited) name("CSV automatic regional encoding `kind'")
}
foreach kind in utf32le utf32be utf16be cp1251 sjis latin2 {
    cio_import_test using "`root'/`kind'.csv", kind(delimited) name("CSV automatic short encoding `kind'")
}
foreach opt in "encoding(bogus)" "delimiters(foo)" "delimiters(|)" {
    capture noisily import delimited using "`root'/empty.csv", `opt' clear
    local nr=_rc
    capture noisily cimport delimited using "`root'/empty.csv", `opt' clear
    local cr=_rc
    cio_record "CSV empty option error `opt'" `=cond(`nr'==`cr',0,9)'
}
di "CSV_OPTIONS_PASSED=$IO_FORMATS_PASSED FAILED=$IO_FORMATS_FAILED"
if $IO_FORMATS_FAILED exit 9
    global CTOOLS_CASE_COMPLETE "io_delimited_options"
    exit
}

if "`__case'" == "io_excel_formulas" {
/* Native/C workbook formula edits: the Python runner also checks XML structure. */
args root
if `"`root'"' == "" local root "`c(pwd)'/temp/io_excel_options"
forvalues test=1/15 {
    local kind simple
    local rows 1
    local cols 1
    local action modify
    local location : word `test' of C2 C3 C4 D2 D3 D4 A2 C2 D3 C2 B2 C2 C2 C2 C2
    if `test'>7 local kind grid
    if `test'==10 | `test'==11 local rows 2
    if `test'==10 local cols 2
    if `test'==11 local cols 3
    if `test'==15 local cols 3
    if `test'==12 {
        local rows 3
        local cols 3
    }
    if `test'==13 local action replace
    if `test'==14 local action new
    clear
    set obs `rows'
    forvalues j=1/`cols' {
        gen double value`j'=99+_n+`j'
    }
    foreach prefix in native result {
        copy "`root'/formulas_`kind'.xlsx" "`root'/formula_`prefix'_`test'.xlsx", replace
        local command export
        if "`prefix'"=="result" local command cexport
        local sheetopt sheet("Data",`action') keepcellfmt
        if "`action'"=="new" local sheetopt sheet("New")
        capture noisily `command' excel using "`root'/formula_`prefix'_`test'.xlsx", `sheetopt' cell(`location')
        local `prefix'rc=_rc
    }
    if `nativerc' | `resultrc' {
        cio_record "XLSX formula edit `test' export" 9
        continue
    }
    local sheets Data
    if "`action'"=="new" local sheets Data New
    foreach sheet of local sheets {
        import excel using "`root'/formula_native_`test'.xlsx", sheet("`sheet'") clear
        local datalabel : data label
        local sortedby : sortedby
        quietly label dir
        local labelsets `r(names)'
        mata: cio_ref=cio_snapshot_data()
        import excel using "`root'/formula_result_`test'.xlsx", sheet("`sheet'") clear
        local datalabel : data label
        local sortedby : sortedby
        quietly label dir
        local labelsets `r(names)'
        mata: cio_compare(cio_ref)
        cio_record "XLSX formula edit `test' `sheet' values/metadata" `=cond(`same',0,9)'
    }
}
    global CTOOLS_CASE_COMPLETE "io_excel_formulas"
    exit
}

if "`__case'" == "io_excel_options" {
/* Differential native options, workbook mutation, and error atomicity. */
do "validation/support/io_test_helpers.do"
args root
if `"`root'"' == "" local root "`c(pwd)'/temp/io_excel_options"
capture mkdir "`root'"
set obs 4
gen double v=_n/7
gen double day=mdy(6,17,2020)+_n
format day %td
gen str6 text="str"
replace text="" in 2
foreach ext in xlsx xls {
    export excel using "`root'/numbers.`ext'", replace firstrow(variables)
}
foreach ext in xls xlsx {
    foreach opt in "allstring" "allstring(%9.2f)" "allstring(%21x)" "allstring(%td)" "allstring(%9.0gc)" "allstring(%9.2e)" "allstring(%9.0g) cellrange(:C4)" "cellrange(B2)" "cellrange(:B4)" "firstrow case(l)" "firstrow case(pre)" {
        cio_import_test using "`root'/numbers.`ext'", kind(excel) name("`ext' `opt'") opts(`opt')
    }
}
cio_import_test using "`root'/numbers", kind(excel) name("Excel extensionless prefers xls")
foreach ext in xlsx xls {
    clear
    set obs 3
    gen double x=_n
    replace x=. in 2
    gen str6 text="str"
    replace text="" in 2
    foreach missing in 99 NA {
        export excel using "`root'/nm.`ext'", missing(`missing') replace
        cexport excel using "`root'/cm.`ext'", missing(`missing') replace
        preserve
        import excel using "`root'/nm.`ext'", clear
        mata: cio_ref=cio_snapshot_data()
        import excel using "`root'/cm.`ext'", clear
        mata: cio_compare(cio_ref)
        local rc=cond(`same',0,9)
        cio_record "`ext' missing(`missing')" `rc'
        restore
    }
}

clear
set obs 4
gen double x=_n
gen str8 s="old"
gen double day=mdy(1,1,2020)+_n
format day %td
foreach ext in xlsx xls {
    export excel using "`root'/base.`ext'", firstrow(variables) sheet(Original) replace
    export excel using "`root'/base.`ext'", firstrow(variables) sheet(Untouched)
}
replace x=10+_n
replace s="new"
replace x=. in 2
replace s="" in 2
foreach ext in xlsx xls {
    foreach opt in "sheet(New)" "sheet(Original)" "sheet(Original, modify)" "sheet(Original, replace)" "sheet(Other, modify)" "sheet(Other, replace)" "sheet(Original, modify) keepcellfmt" "sheet(Original) sheetmodify" "sheet(Original) sheetreplace" "sheet(Original, modify) replace" "keepcellfmt" "sheet(Original, replace) keepcellfmt" "sheet(,modify)" "sheet(Original, modify) missing(99)" "sheet(Original, modify) missing(NA)" "sheet(Original, modify) firstrow(var)" {
        copy "`root'/base.`ext'" "`root'/native.`ext'", replace
        copy "`root'/base.`ext'" "`root'/result.`ext'", replace
        capture noisily export excel using "`root'/native.`ext'", `opt' cell(B2)
        local nr=_rc
        capture noisily cexport excel using "`root'/result.`ext'", `opt' cell(B2)
        local cr=_rc
        if `nr'!=`cr' cio_record "`ext' edit error `opt'" 9
        else if `nr' cio_record "`ext' edit error `opt'" 0
        else {
            preserve
            import excel using "`root'/native.`ext'", describe
            local k=r(N_worksheet)
            forvalues i=1/`k' {
                local sh`i' `"`r(worksheet_`i')'"'
            }
            forvalues i=1/`k' {
                import excel using "`root'/native.`ext'", sheet(`"`sh`i''"') clear
                local datalabel : data label
                local sortedby : sortedby
                quietly label dir
                local labelsets `r(names)'
                mata: cio_ref=cio_snapshot_data()
                import excel using "`root'/result.`ext'", sheet(`"`sh`i''"') clear
                local datalabel : data label
                local sortedby : sortedby
                quietly label dir
                local labelsets `r(names)'
                mata: cio_compare(cio_ref)
                local rc=cond(`same',0,9)
                cio_record "`ext' edit `opt' `sh`i''" `rc'
            }
            restore
        }
    }
}
clear
set obs 2
forvalues j=1/8 {
    gen v`j'=_n
}
label var v1 "A"
label var v2 "A"
label var v3 "123"
label var v4 "if"
label var v5 "A B"
label var v6 "C"
label var v7 "ΩABC"
label var v8 "ΩABC"
foreach ext in xls xlsx {
    export excel using "`root'/headers.`ext'", firstrow(varlabels) replace
}
foreach ext in xls xlsx {
    cio_import_test using "`root'/headers.`ext'", kind(excel) name("`ext' duplicate/numeric/reserved/Unicode names") opts(firstrow)
}
foreach opt in "allstring allstring(%9.2f)" "allstring(%9s)" "allstring(bad)" "case(p)" "cellrange(B2:)" {
    capture noisily import excel using "`root'/numbers.xlsx", `opt' clear
    local nr=_rc
    capture noisily cimport excel using "`root'/numbers.xlsx", `opt' clear
    local cr=_rc
    local rc=cond(`nr'==`cr',0,9)
    cio_record "invalid Excel import option `opt'" `rc'
}
clear
set obs 3
gen byte x=_n
cio_export_test, kind(dbase) ext(dbf) name("DBF version(3)") opts(version(3))
foreach ver in 4 III IV iii iv {
    capture noisily export dbase using "`root'/native.dbf", version(`ver') replace
    local nr=_rc
    capture noisily cexport dbase using "`root'/result.dbf", version(`ver') replace
    local cr=_rc
    local rc=cond(`nr'==`cr',0,9)
    cio_record "DBF rejects version(`ver')" `rc'
}
foreach ext in xls xlsx {
    foreach sel in "id=A note=C" "id=B note=D" "empty=D" "x=A y=A" "x y" "x y z w" {
        foreach opt in "" "allstring" "allstring(%9.2f)" {
            cio_import_test using "`root'/numbers.`ext'", kind(excel) name("`ext' selector `sel' `opt'") select(`sel') opts(`opt')
        }
    }
    foreach opt in "locale(UTF-8)" "locale(ISO-8859-1)" "locale(en_US)" "locale(tr_TR)" "locale(bogus)" {
        cio_import_test using "`root'/headers.`ext'", kind(excel) name("`ext' Unicode header `opt'") opts(firstrow `opt')
    }
    capture noisily import excel id=A using "`root'/numbers.`ext'", firstrow clear
    local nr=_rc
    capture noisily cimport excel id=A using "`root'/numbers.`ext'", firstrow clear
    local cr=_rc
    cio_record "`ext' extvarlist conflicts with firstrow" `=cond(`nr'==`cr',0,9)'
    foreach sel in "empty=Z" "id=3" "x=A x=B" {
        capture noisily import excel `sel' using "`root'/numbers.`ext'", clear
        local nr=_rc
        capture noisily cimport excel `sel' using "`root'/numbers.`ext'", clear
        local cr=_rc
        cio_record "`ext' invalid selector `sel'" `=cond(`nr'==`cr',0,9)'
    }
}
foreach ext in xls xlsx {
    foreach opt in "firstrow" "firstrow allstring" "firstrow allstring(%9.2f)" {
        cio_import_test using "`root'/display_formats.`ext'", kind(excel) name("`ext' built-in/custom date and time formats `opt'") opts(`opt')
    }
}
do "validation/stata/io.do" io_excel_formulas "`root'"
foreach ext in xls xlsx {
    foreach opt in "clear" "firstrow" "allstring" "sheet(Formats)" "cellrange(A1)" "detail" {
        capture noisily import excel using "`root'/display_formats.`ext'", describe `opt'
        local nr=_rc
        capture noisily cimport excel using "`root'/display_formats.`ext'", describe `opt'
        cio_record "`ext' describe rejects `opt'" `=cond(`nr'==_rc,0,9)'
    }
    foreach corruption in invalid truncated {
        cio_import_state_test using "`root'/`ext'_`corruption'.`ext'", kind(excel) state(changed) opts(describe)
    }
}
di "IO_OPTIONS_PASSED=$IO_FORMATS_PASSED FAILED=$IO_FORMATS_FAILED"
if $IO_FORMATS_FAILED exit 9
    global CTOOLS_CASE_COMPLETE "io_excel_options"
    exit
}

if "`__case'" == "io_formats" {
do "validation/support/io_test_helpers.do" `0'

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
foreach year in 1987 2000 2025 {
    // DBF metadata stores years since 1900, as native import dbase does.
    char _dta[dbf_year] `=`year'-1900'
    char _dta[dbf_month] 2
    char _dta[dbf_day] 17
    cio_export_test, kind(dbase) ext(dbf) name("dBase original header date `year'") opts(origdbfdate)
    mata: fh=fopen("temp/io_allformats/native/output.dbf","r"); ref_header=fread(fh,4); fclose(fh)
    mata: fh=fopen("temp/io_allformats/result/output.dbf","r"); got_header=fread(fh,4); fclose(fh)
    mata: st_local("header_same",strofreal(ref_header==got_header))
    cio_record "dBase original date bytes `year'" `=cond(`header_same',0,9)'
}

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
    else {
        cio_export_test, kind(shp) ext(shp) name("SHP `shape' writer") opts(shx)
        clonevar custom_id=_ID
        replace custom_id=custom_id+100 if !missing(custom_id)
        cio_export_test, kind(shp) ext(shp) name("SHP `shape' explicit identifier") opts(shx id(custom_id))
    }
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
    global CTOOLS_CASE_COMPLETE "io_formats"
    exit
}

if "`__case'" == "io_parity" {
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
    global CTOOLS_CASE_COMPLETE "io_parity"
    exit
}

di as error "Unknown io case: `__case'"
exit 198
