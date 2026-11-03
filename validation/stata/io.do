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
* Lone CR line endings end rows as in native: classic Mac files, mixed CR and
* LF, CR inside quotes, blank CR lines and CR before CRLF.
foreach pair in `"cr_only id,value<CR>1,100<CR>2,200<CR>"' `"cr_mixed id,value<CR>1,100<LF>2,200<CR>3,300<LF>"' `"cr_quoted id,text<CR>1,"a<CR>b"<CR>2,"c"<CR>"' `"cr_blank id,value<CR><CR>1,100<CR><CR>2,200<CR>"' `"cr_crlf id,value<CR><CR><LF>1,100<CR><CR><LF>2,200<CR><CR><LF>"' `"cr_noend id<TAB>value<CR>1<TAB>100<CR>2<TAB>200"' {
    gettoken file text : pair
    mata: _p = st_local("root") + "/`file'.csv"; _t = subinstr(subinstr(subinstr(strtrim(st_local("text")), "<CR>", char(13)), "<LF>", char(10)), "<TAB>", char(9)); if (fileexists(_p)) unlink(_p); _f = fopen(_p, "w"); fwrite(_f, _t); fclose(_f)
    foreach opt in "" "bindquotes(strict)" "rowrange(2:3)" "varnames(nonames)" "bindquotes(nobind)" {
        cio_import_test using "`root'/`file'.csv", kind(delimited) name("CSV lone CR `file' `opt'") opts(`opt')
    }
}
* colrange() selects columns before header detection, naming and
* numericcols()/stringcols() numbering, as in native.
foreach pair in "cr_dup id,x,x<LF>1,2,3<LF>4,5,6<LF>" "cr_hdr a,b<LF>1,x<LF>2,y<LF>" "cr_force a,b,c<LF>1,2,3<LF>4,x,6<LF>" "cr_force2 a,b,c<LF>q,2,3<LF>4,5,6<LF>" "cr_ragged a,b<LF>1,2,3<LF>4,5<LF>" "cr_nohdr 1,2,3<LF>4,5,6<LF>" "cr_emptyhdr a,,c<LF>1,2,3<LF>" {
    gettoken file text : pair
    mata: _p = st_local("root") + "/`file'.csv"; _t = subinstr(subinstr(subinstr(strtrim(st_local("text")), "<CR>", char(13)), "<LF>", char(10)), "<TAB>", char(9)); if (fileexists(_p)) unlink(_p); _f = fopen(_p, "w"); fwrite(_f, _t); fclose(_f)
}
foreach case in "dup colrange(3)" "dup colrange(2:3)" "dup colrange(2:3) varnames(1)" "hdr colrange(2)" "hdr colrange(2) varnames(1)" "hdr colrange(3:4)" "hdr colrange(5)" "force colrange(2:3) numericcols(1)" "force colrange(2:3) stringcols(3)" "force2 colrange(2:3) stringcols(1)" "ragged colrange(2)" "ragged colrange(3)" "nohdr colrange(2:3)" "nohdr colrange(2:3) varnames(nonames)" "emptyhdr colrange(2:3)" {
    gettoken file opt : case
    cio_import_test using "`root'/cr_`file'.csv", kind(delimited) name("CSV colrange `file' `opt'") opts(`opt')
}
* A file wider than maxvar imports with a narrow colrange().
local file cr_wide
mata: st_local("text", invtokens(strofreal(1..(c("maxvar") + 10)), ",") + "<LF>" + invtokens(strofreal(1..(c("maxvar") + 10)), ",") + "<LF>")
mata: _p = st_local("root") + "/`file'.csv"; _t = subinstr(subinstr(subinstr(strtrim(st_local("text")), "<CR>", char(13)), "<LF>", char(10)), "<TAB>", char(9)); if (fileexists(_p)) unlink(_p); _f = fopen(_p, "w"); fwrite(_f, _t); fclose(_f)
cio_import_test using "`root'/cr_wide.csv", kind(delimited) name("CSV colrange on a file wider than maxvar") opts(colrange(1:10) varnames(nonames))
mata: mata drop _p _t _f
* Quotes as native reads them: only a quote at a field's first byte opens it
* and is dropped (with a closing quote in its last byte); "" is one quote
* inside quotes only; without binding a run of quotes is one; under
* bindquotes(strict) a CRLF inside quotes keeps only the CR.
foreach opt in "" "bindquotes(nobind)" "stripquotes(yes)" "stripquotes(no)" {
    cio_import_test using "`root'/quote_forms.csv", kind(delimited) name("CSV quote forms `opt'") opts(stringcols(2) `opt')
}
foreach opt in "" "bindquotes(strict)" "bindquotes(nobind)" {
    cio_import_test using "`root'/quote_numbers_padded.csv", kind(delimited) name("CSV space-padded quoted numbers `opt'") opts(`opt')
    cio_import_test using "`root'/quote_runs.csv", kind(delimited) name("CSV doubled and run quotes `opt'") opts(`opt')
}
foreach opt in "bindquotes(strict)" "bindquotes(strict) stripquotes(yes)" "bindquotes(strict) stripquotes(no)" {
    cio_import_test using "`root'/quote_crlf.csv", kind(delimited) name("CSV CRLF inside quotes `opt'") opts(`opt')
}
* NaN and Infinity (exact case, signed or not) are missing values, and d/D
* exponents read as e; nan, NA and malformed exponents stay text.
cio_import_test using "`root'/number_tokens.csv", kind(delimited) name("CSV NaN, Infinity and d exponents")
cio_import_test using "`root'/number_strings.csv", kind(delimited) name("CSV nan, NA and malformed exponents stay text")
* numericcols() wins over stringcols(), including _all.
foreach opt in "stringcols(2 3) numericcols(3)" "stringcols(_all) numericcols(2)" "numericcols(_all) stringcols(2)" "stringcols(2) numericcols(_all)" {
    cio_import_test using "`root'/force.csv", kind(delimited) name("CSV overlapping overrides `opt'") opts(`opt')
}
* rowrange(): variables come from the first row (unless nonames), the
* selected rows and the row that seeds an explicit header's inference.
foreach file in a b c d {
    foreach opt in "rowrange(3:4)" "rowrange(2:2)" "rowrange(4)" "rowrange(3:4) varnames(1)" "rowrange(2:3) varnames(nonames)" "rowrange(1:2) varnames(nonames)" "rowrange(:2)" "rowrange(10:12)" "rowrange(2:3)" "rowrange(3:3) varnames(nonames)" {
        cio_import_test using "`root'/ragged_`file'.csv", kind(delimited) name("CSV ragged rows under `opt' (`file')") opts(`opt')
    }
}
cio_import_test using "`root'/long_ids.csv", kind(delimited) name("CSV long variables display as %12.0g")
* Display formats follow the type native infers from its first 513 lines.
foreach tag in long int double float {
    foreach row in 100 512 513 5000 {
        cio_import_test using "`root'/sample_`tag'_`row'.csv", kind(delimited) name("CSV format from the first 513 lines, `tag' outlier at row `row'")
    }
}
* Header names: Unicode letters, decimal digits and underscores, no leading
* digits, 32 characters, then the case (a case-changed name over 32
* characters falls back to v#/V#).
foreach case in preserve lower upper {
    cio_import_test using "`root'/unicode_headers.csv", kind(delimited) name("CSV Unicode header names `case'") opts(case(`case') encoding(utf8))
}
cio_import_test using "`root'/wide_names.csv", kind(delimited) name("CSV 2,000 columns with 32-character names")
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
* Header fallbacks can displace earlier names, including collision chains.
foreach headers in "B B C" "B C C" "B C B" "C C A" "B 123 C" "B B A" "X X X" "A A A" {
    clear
    set obs 2
    forvalues j=1/3 {
        gen double v`j'=`j'+_n/7
        local header : word `j' of `headers'
        label var v`j' "`header'"
    }
    foreach ext in xls xlsx {
        export excel using "`root'/collisions.`ext'", firstrow(varlabels) replace
    }
    foreach ext in xls xlsx {
        foreach case in preserve lower upper {
            cio_import_test using "`root'/collisions.`ext'", kind(excel) name("`ext' fallback chain `headers' `case'") opts(firstrow case(`case'))
        }
    }
}
* Native resolves raw names before case conversion; a new collision is an
* error. Our error contract also preserves the caller's original dataset.
clear
set obs 2
gen a=_n
gen b=_n+1
label var a "x"
label var b "X"
export excel using "`root'/case_collision.xlsx", firstrow(varlabels) replace
capture noisily import excel using "`root'/case_collision.xlsx", firstrow case(lower) clear
local native_rc=_rc
clear
set obs 2
gen double sentinel=99
label var sentinel "Original label"
format sentinel %12.0g
capture noisily cimport excel using "`root'/case_collision.xlsx", firstrow case(lower) clear
assert _rc==`native_rc' & `native_rc'==110
assert _N==2 & c(k)==1
assert sentinel==99
assert "`: variable label sentinel'"=="Original label"
assert "`: format sentinel'"=="%12.0g"
cio_record "Excel case collision preserves original data" 0
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
* cellrange() must lie within the worksheet's data area, as in native. Column
* B has a value only below the range, so A1:B2 is still inside it.
clear
set obs 5
gen a=_n
gen b=cond(_n==5,9,.)
quietly export excel using "`root'/late_column.xls", replace
quietly export excel using "`root'/late_column.xlsx", replace
foreach ext in xls xlsx {
    cio_import_test using "`root'/late_column.`ext'", kind(excel) name("`ext' cellrange inside a later data column") opts(cellrange(A1:B2))
    foreach range in "A1:C2" "A1:B6" "C1:C2" "A6:B6" "A1:IV3" "A1:A65536" "C3" {
        capture noisily import excel using "`root'/late_column.`ext'", clear cellrange(`range')
        local nr=_rc
        capture noisily cimport excel using "`root'/late_column.`ext'", clear cellrange(`range')
        cio_record "`ext' cellrange(`range') outside the data area" `=cond(`nr'==_rc & _rc==198,0,9)'
    }
}
* More columns than maxvar, and a lone value at XFD1048576 (17 billion dense
* cells), fail before the data in memory is replaced.
quietly putexcel set "`root'/wide.xlsx", replace
quietly putexcel A1 = 1
quietly putexcel XFD1 = 2
quietly putexcel save
quietly putexcel set "`root'/far.xlsx", replace
quietly putexcel A1 = 1
quietly putexcel XFD1048576 = 2
quietly putexcel save
capture noisily import excel using "`root'/wide.xlsx", clear
local nr=_rc
foreach book in wide far {
    clear
    quietly set obs 2
    gen long sentinel=99
    capture noisily cimport excel using "`root'/`book'.xlsx", clear
    local cr=_rc
    local kept=c(k)==1 & _N==2 & sentinel[1]==99
    if "`book'"=="wide" cio_record "xlsx columns beyond maxvar keep the data" `=cond(`nr'==900 & `cr'==900 & `kept',0,9)'
    else cio_record "xlsx cell at XFD1048576 fails with r(909)" `=cond(`cr'==909 & `kept',0,9)'
}
* sheet(..., modify) on workbooks whose rows or cells omit their optional r
* attributes, and on one with no styles part: compare with native's edit.
foreach book in norcells norboth nostyles {
    foreach who in export cexport {
        copy "`root'/positions_`book'.xlsx" "`root'/positions_`book'_`who'.xlsx", replace
        clear
        quietly set obs 3
        gen x=9000+_n
        quietly `who' excel x using "`root'/positions_`book'_`who'.xlsx", sheet("Sheet1", modify) cell(B2)
    }
    quietly import excel using "`root'/positions_`book'_export.xlsx", clear
    tempfile edited
    quietly save `edited'
    capture noisily import excel using "`root'/positions_`book'_cexport.xlsx", clear
    local rc=_rc
    if !`rc' {
        capture cf _all using `edited', all
        local rc=_rc
    }
    cio_record "xlsx modify of positions_`book' matches native" `rc'
}
* Text that is not valid UTF-8 fails like native, keeping the destination.
clear
quietly set obs 2
gen str10 s=cond(_n==1,"caf"+char(233),"ok")
capture noisily export excel using "`root'/invalid_utf8_native.xlsx", replace
local nr=_rc
quietly export excel using "`root'/invalid_utf8.xlsx" in 2, replace
quietly checksum "`root'/invalid_utf8.xlsx"
local before=r(checksum)
foreach mode in replace modify {
    local opt replace
    if "`mode'"=="modify" local opt `"sheet("Sheet1", modify)"'
    capture noisily cexport excel using "`root'/invalid_utf8.xlsx", `opt'
    local cr=_rc
    quietly checksum "`root'/invalid_utf8.xlsx"
    cio_record "xlsx invalid UTF-8 `mode' fails like native" `=cond(`nr'==9901 & `cr'==9901 & r(checksum)==`before',0,9)'
}
* XLS modify of a worksheet that embeds a chart substream: the chart is kept
* and both sheets stay readable. (Native's own edit of this file empties them.)
copy "`root'/chart_in_sheet.xls" "`root'/chart_in_sheet_edit.xls", replace
clear
quietly set obs 2
gen c=_n+6
capture noisily cexport excel c using "`root'/chart_in_sheet_edit.xls", sheet("Data", modify) cell(D2)
local rc=_rc
if !`rc' {
    capture noisily import excel using "`root'/chart_in_sheet_edit.xls", sheet("Data") clear
    local rc=_rc
    if !`rc' local rc=cond(_N==4 & c(k)==4 & A[2]==11 & B[4]==32 & D[2]==7 & D[3]==8 & missing(D[1]),0,9)
}
if !`rc' {
    capture noisily import excel using "`root'/chart_in_sheet_edit.xls", sheet("Other") clear
    local rc=cond(_rc,_rc,cond(_N==3 & A[3]==102,0,9))
}
cio_record "xls modify keeps an embedded chart substream" `rc'
* A chart sheet first: cimport reads the first worksheet, rejects sheet("Chart1")
* (error 601) for import and export, and describe lists it without a range.
quietly import excel using "`root'/chart_sheet_first.xls", sheet("Data") clear
tempfile first_worksheet
quietly save `first_worksheet'
capture noisily cimport excel using "`root'/chart_sheet_first.xls", clear
local rc=_rc
if !`rc' {
    capture cf _all using `first_worksheet', all
    local rc=_rc
}
cio_record "xls default sheet skips a chart sheet" `rc'
capture noisily cimport excel using "`root'/chart_sheet_first.xls", sheet("Chart1") clear
cio_record "xls sheet() naming a chart sheet fails" `=cond(_rc==601,0,9)'
cio_workbook_test using "`root'/chart_sheet_first.xls", name("xls describe lists a chart sheet without a range")
quietly checksum "`root'/chart_sheet_first.xls"
local before=r(checksum)
clear
quietly set obs 1
gen c=9
capture noisily cexport excel c using "`root'/chart_sheet_first.xls", sheet("Chart1", modify) cell(B2)
local cr=_rc
quietly checksum "`root'/chart_sheet_first.xls"
cio_record "xls modify of a chart sheet fails and keeps the file" `=cond(`cr'==601 & r(checksum)==`before',0,9)'
* Header names as native: Unicode letters, decimal digits and underscores of
* the first 32 characters, leading digits dropped, validity judged in the
* requested case (so _cons and str10 stand under upper).
clear
quietly set obs 2
local headers `"1st_name|2019_revenue|3M|9 lives|12|٣x|a|1a|Ünïcödé_ñame|école|Straße|ǅb|µm|x²y|a€b|mad½z|ßßßßßßßßßßßßßßßß|1234567890bcdefghijklmnopqrstuvwxyzABC|€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€xyz"'
mata: _h = tokens(st_local("headers"), "|"); _h = select(_h, _h :!= "|")
mata: for (_j = 1; _j <= cols(_h); _j++) (void) st_addvar("double", "v" + strofreal(_j))
mata: for (_j = 1; _j <= cols(_h); _j++) st_varlabel(_j, _h[_j])
mata: mata drop _h _j
foreach v of varlist _all {
    quietly replace `v' = _n
}
foreach ext in xls xlsx {
    export excel using "`root'/unicode_headers.`ext'", firstrow(varlabels) replace
    foreach case in preserve lower upper {
        cio_import_test using "`root'/unicode_headers.`ext'", kind(excel) name("`ext' Unicode and digit header names `case'") opts(firstrow case(`case'))
    }
}
* Native creates _cons and str# headers as named (an error, r(198)), but they
* stand as _CONS and STR10 under case(upper).
foreach header in _cons str10 {
    foreach ext in xls xlsx {
        clear
        quietly set obs 2
        gen double v1 = _n
        gen double v2 = _n
        label var v1 "`header'"
        label var v2 "other"
        quietly export excel using "`root'/reserved_header.`ext'", firstrow(varlabels) replace
        cio_import_test using "`root'/reserved_header.`ext'", kind(excel) name("`ext' header `header' under case(upper)") opts(firstrow case(upper))
        foreach case in preserve lower {
            capture noisily import excel using "`root'/reserved_header.`ext'", firstrow case(`case') clear
            local nr=_rc
            clear
            quietly set obs 1
            gen byte sentinel=7
            capture noisily cimport excel using "`root'/reserved_header.`ext'", firstrow case(`case') clear
            local cr=_rc
            capture assert c(k)==1 & sentinel==7
            cio_record "`ext' header `header' fails as native under case(`case')" `=cond(`nr'==198 & `cr'==198 & !_rc,0,9)'
        }
    }
}
* Long sheet names, entity-heavy names, more than 256 sheets, and cells without
* a value at the edges of the used range (excel_import_fixtures).
local cjk "東日本営業本部売上実績四半期報告書第三版改訂済み最終"
mata: st_local("amp", invtokens(tokens("A B C D E F G H I J K L M N O P"), "&"))
foreach s in `"`cjk'"' `"`amp'"' "Plain" {
    cio_import_test using "`root'/sheet_names.xlsx", kind(excel) name("xlsx sheet() with a long or entity-heavy name") opts(sheet(`"`s'"') firstrow)
}
cio_workbook_test using "`root'/sheet_names.xlsx", name("xlsx describe long and entity-heavy sheet names")
foreach s in Sheet1 Sheet256 Sheet257 Sheet300 {
    cio_import_test using "`root'/many_sheets.xlsx", kind(excel) name("xlsx sheet(`s') of 300 sheets") opts(sheet(`s'))
}
cio_workbook_test using "`root'/many_sheets.xlsx", name("xlsx describe 300 sheets")
foreach f in styled rowcells emptyv emptystr emptyis formula error errall leadcol {
    foreach opt in "" "firstrow" "allstring" {
        cio_import_test using "`root'/gap_`f'.xlsx", kind(excel) name("xlsx used range with valueless cells `f' `opt'") opts(`opt')
    }
    cio_workbook_test using "`root'/gap_`f'.xlsx", name("xlsx describe used range with valueless cells `f'")
}
* Malformed and reversed cellrange() corners fail with r(198), keeping the data.
foreach ext in xls xlsx {
    foreach range in "A1:C" "A1:D2000000" "1:10" "D10:A1" "A0:C5" "XFE1:XFF3" "A1:C5x" "B2:A3" {
        capture noisily import excel using "`root'/late_column.`ext'", clear cellrange(`range')
        local nr=_rc
        clear
        quietly set obs 1
        gen byte sentinel=7
        capture noisily cimport excel using "`root'/late_column.`ext'", clear cellrange(`range')
        local cr=_rc
        capture assert c(k)==1 & sentinel==7
        cio_record "`ext' malformed cellrange(`range') fails and keeps the data" `=cond(`nr'==198 & `cr'==198 & !_rc,0,9)'
    }
}
* ZIP part sizes are checked before anything is allocated: a declared size
* its deflate data cannot reach fails as native (r(603)), and a part beyond
* c(max_memory) fails with r(909). Neither replaces the data in memory.
capture noisily import excel using "`root'/zip_declared.xlsx", clear
local nr=_rc
clear
quietly set obs 1
gen byte sentinel=7
capture noisily cimport excel using "`root'/zip_declared.xlsx", clear
local cr=_rc
capture assert c(k)==1 & sentinel==7
cio_record "xlsx part declaring a false size fails as native" `=cond(`nr'==603 & `cr'==603 & !_rc,0,9)'
clear
quietly set max_memory 100m
quietly set obs 1
gen byte sentinel=7
capture noisily cimport excel using "`root'/zip_150mb.xlsx", clear
local cr=_rc
capture assert c(k)==1 & sentinel==7
local kept=!_rc
clear
quietly set max_memory .
capture noisily cimport excel using "`root'/zip_150mb.xlsx", clear
cio_record "xlsx part beyond max_memory fails with r(909)" `=cond(`cr'==909 & `kept' & !_rc & _N==3,0,9)'
* A chart embedded in an XLS worksheet keeps its cached values to itself.
foreach sheet in Data Other {
    cio_import_test using "`root'/chart_in_sheet.xls", kind(excel) name("xls worksheet with an embedded chart `sheet'") opts(sheet(`sheet'))
}
cio_workbook_test using "`root'/chart_in_sheet.xls", name("xls describe a worksheet with an embedded chart")
* Numbers in string columns read as native's %18.0g text; numeric header
* cells label their variables with that text.
foreach ext in xls xlsx {
    foreach opt in "" "allstring" "firstrow" "firstrow allstring" {
        cio_import_test using "`root'/string_numbers.`ext'", kind(excel) name("`ext' numbers in a string column `opt'") opts(`opt')
    }
}
* Value text: a non-number reads as 0 and a number prefix as that number;
* a formula text result decodes _xHHHH_ escapes.
foreach opt in "" "allstring" {
    cio_import_test using "`root'/value_text.xlsx", kind(excel) name("xlsx cell value text `opt'") opts(`opt')
}
cio_import_test using "`root'/numeric_headers.xlsx", kind(excel) name("xlsx numeric header labels") opts(firstrow)
/* cexport excel options compared with export excel: rc and every worksheet
   (name, range and contents) of the resulting workbook. */
capture program drop _g4b_book
program define _g4b_book
    args file
    mata: _g4b_state = J(1, 0, "")
    capture quietly import excel using "`file'", describe
    if _rc exit
    local n = r(N_worksheet)
    forvalues j = 1/`n' {
        local s`j' `"`r(worksheet_`j')'"'
        mata: _g4b_state = _g4b_state, st_global("r(worksheet_`j')"), st_global("r(range_`j')")
    }
    forvalues j = 1/`n' {
        capture quietly import excel using "`file'", sheet(`"`s`j''"') allstring clear
        if !_rc & _N {
            mata: _g4b_state = _g4b_state, invtokens(vec(st_sdata(., .))')
        }
    }
end
capture program drop _g4b_excel
program define _g4b_excel
    syntax, root(string) name(string) [opts(string asis) setup(string asis) ext(string) vars(string) sel(string asis)]
    if "`ext'" == "" local ext xlsx
    * cio_record displays the name in double quotes
    mata: st_local("name", subinstr(st_local("name"), char(34), char(39)))
    tempfile data
    quietly save `data'
    foreach side in native result {
        local cmd = cond("`side'" == "native", "export", "cexport")
        local f "`root'/g4b_`side'.`ext'"
        capture erase "`f'"
        if `"`setup'"' != "" {
            quietly export excel using "`f'", sheet("a") replace
            quietly export excel using "`f'", sheet("Sheet1")
            quietly keep in 1/2
            quietly replace price = -1
        }
        capture noisily `cmd' excel `vars' using "`f'" `sel', `opts'
        local rc_`side' = _rc
        _g4b_book "`f'"
        mata: _g4b_`side' = _g4b_state
        use `data', clear
    }
    mata: st_local("same", strofreal(_g4b_native == _g4b_result))
    cio_record "`name' (native rc=`rc_native', C rc=`rc_result')" `=cond(`rc_native' == `rc_result' & `same', 0, 9)'
end
sysuse auto, clear
quietly keep in 1/3
keep make price
foreach o in `"sheet("a", modify) sheetreplace"' `"sheet("a", replace) sheetmodify"' `"sheet("a", modify replace)"' `"sheet("a", replace modify)"' `"sheet("a", rep mod)"' `"sheet("a", modify modify)"' "sheetmodify" "sheetreplace" "sheetmodify sheetreplace" "sheetreplace replace" {
    _g4b_excel, root("`root'") name(`"excel `o'"') opts(`o') setup(yes)
}
foreach o in `"sheet(`"a" b"')"' `"sheet(`""a""')"' `"sheet(a b)"' `"sheet("a" b)"' `"sheet("a b")"' `"sheet(a, bogus)"' {
    _g4b_excel, root("`root'") name(`"excel `o'"') opts(`o')
}
/* A value-labeled variable with a 32-character name */
clear
quietly set obs 3
gen byte respondent_employment_status_cat = _n
label define g4b_emp 1 "employed" 2 "unemployed"
label values respondent_employment_status_cat g4b_emp
_g4b_excel, root("`root'") name("excel labeled 32-character name") opts(replace firstrow(variables))
/* Zero observations: r(2000); an empty selection: r(198) */
sysuse auto, clear
foreach ext in xlsx xls {
    _g4b_excel, root("`root'") name("`ext' empty selection") vars(price) sel(if price < 0) opts(replace) ext(`ext')
    preserve
    quietly drop in 1/l
    _g4b_excel, root("`root'") name("`ext' zero observations") vars(price) opts(replace) ext(`ext')
    _g4b_excel, root("`root'") name("`ext' zero observations firstrow") vars(price) opts(replace firstrow(variables)) ext(`ext')
    restore
}
/* An existing file that is not a workbook, without replace or with modify */
foreach ext in xlsx xls {
    foreach o in "" `"sheet("x", modify)"' {
        local what = cond(`"`o'"' == "", "", "with sheet(x, modify)")
        foreach side in native result {
            local cmd = cond("`side'" == "native", "export", "cexport")
            capture erase "`root'/notbook_`side'.`ext'"
            file open fh using "`root'/notbook_`side'.`ext'", write text replace
            file write fh "not a workbook" _n
            file close fh
            capture noisily `cmd' excel price using "`root'/notbook_`side'.`ext'", `o'
            local rc_`side' = _rc
        }
        cio_record "`ext' existing non-workbook `what' (native rc=`rc_native', C rc=`rc_result')" `=cond(`rc_native' == `rc_result', 0, 9)'
    }
}
/* More than 256 columns to .xls */
clear
quietly set obs 2
forvalues i = 1/257 {
    quietly gen byte v`i' = `i'
}
_g4b_excel, root("`root'") name("xls 257 columns") opts(replace) ext(xls)
/* keepcellfmt: numbers keep the cell's format, dates get a date format */
foreach ext in xlsx xls {
    foreach side in native result {
        local cmd = cond("`side'" == "native", "export", "cexport")
        clear
        quietly set obs 5
        forvalues j = 1/4 {
            quietly gen double c`j' = _n * `j' + 0.5
        }
        quietly export excel using "`root'/keep_`side'.`ext'", replace
        clear
        quietly set obs 3
        gen double n = _n + 0.25
        gen d = td(01mar2021) + _n
        format d %td
        gen double t = clock("01jan2021 01:02:03", "DMYhms") + _n
        format t %tc
        capture noisily `cmd' excel using "`root'/keep_`side'.`ext'", sheet("Sheet1", modify) cell(B2) keepcellfmt
        local rc_`side' = _rc
        quietly import excel using "`root'/keep_`side'.`ext'", clear
        mata: _g4b_`side' = J(1, 0, "")
        foreach v of varlist _all {
            mata: _g4b_`side' = _g4b_`side', st_varformat("`v'"), invtokens(strofreal(st_data(., "`v'"), "%21x")')
        }
    }
    mata: st_local("same", strofreal(_g4b_native == _g4b_result))
    cio_record "`ext' keepcellfmt with dates (native rc=`rc_native', C rc=`rc_result')" `=cond(`rc_native' == `rc_result' & `same', 0, 9)'
}
/* %tb converts through its business calendar; %tg is written as its number.
   (Native writes %tb as empty text and %tg as a datetime.) */
file open fh using "`root'/g4bcal.stbcal", write text replace
file write fh "version 14" _n "purpose " _char(34) "Mon-Fri 2018" _char(34) _n "dateformat dmy" _n ///
    "range 01jan2018 31dec2018" _n "centerdate 01jan2018" _n "omit dayofweek (Sa Su)" _n
file close fh
adopath ++ "`root'"
foreach ext in xlsx xls {
    clear
    quietly set obs 3
    gen double day = mdy(1, 1, 2018) + 7 * (_n - 1)
    gen double bd = bofd("g4bcal", day)
    format bd %tbg4bcal
    gen double bd2 = bd
    format bd2 %tbg4bcal:CCYY-NN-DD
    gen double g = _n + 0.5
    format g %tg
    capture noisily cexport excel bd bd2 g using "`root'/tb.`ext'", replace firstrow(variables)
    local rc = _rc
    if !`rc' {
        tempfile orig
        quietly save `orig'
        quietly import excel using "`root'/tb.`ext'", firstrow clear
        quietly merge 1:1 _n using `orig', keepusing(day) nogenerate
        local fbd : format bd
        local fg : format g
        quietly count if bd != day | bd2 != day | g != _n + 0.5
        local rc = cond(r(N) == 0 & substr("`fbd'", 1, 3) == "%td" & substr("`fg'", 1, 3) != "%td" & substr("`fg'", 1, 3) != "%tc", 0, 9)
    }
    cio_record "`ext' %tb through its calendar, %tg as a number" `rc'
}
adopath - "`root'"
capture mata: mata drop _g4b_state _g4b_native _g4b_result
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
preserve
clear
set obs 6
gen double day=cond(_n==1,td(01jan1960),cond(_n==2,td(01jan1973),td(01jan2020)))
gen double leftday=day
gen double clock=day*86400000+43200000+.125
gen double leftclk=clock
gen double leap=Cofc(clock)
gen double leftleap=leap
foreach v of varlist _all {
    replace `v'=. in 5
    replace `v'=.z in 6
    label var `v' "Preserved `v' label"
}
format day %td
format leftday %-td
format clock %tc
format leftclk %-tc
format leap %tC
format leftleap %-tC
cio_export_test, kind(spss) ext(sav) name("SPSS left-justified dates and leap-second clocks")
cio_export_test, kind(spss) ext(sav) name("SPSS selected leap-second clocks") select(leap leftleap if _n!=2)
* The source values and formats must survive both success and publication failure.
assert leap==Cofc(clock) if !missing(clock)
assert "`: format leftleap'"=="%-tC"
capture noisily cexport spss using "temp/io_allformats/result/output.sav", replace
assert _rc==0
capture noisily cexport spss using "temp/io_allformats/result/output.sav"
assert _rc==602
assert leap==Cofc(clock) if !missing(clock)
assert leap==.z & leftleap==.z in 6
foreach pair in "sasxport5 xpt" "sasxport8 v8xpt" {
    gettoken kind ext : pair
    local ext=strtrim("`ext'")
    cio_export_test, kind(`kind') ext(`ext') name("`kind' left date formats retain transport units") select(day leftday clock leftclk)
}
* Native export writes the legacy %d format as a plain number, not a date.
gen double legacy=day
gen double leftlegacy=day
format legacy %d
format leftlegacy %-d
foreach pair in "spss sav" "sasxport8 v8xpt" "dbase dbf" {
    gettoken kind ext : pair
    local ext=strtrim("`ext'")
    cio_export_test, kind(`kind') ext(`ext') name("`kind' legacy %d format is numeric") select(legacy leftlegacy)
}
* Native SPSS and XPORT8 writers store period dates as the first day of the
* period with a period format, and store %tC clocks without leap seconds.
gen double month=cond(missing(day),day,mofd(day)+_n)
replace month=month+.5 in 4
gen double leftmonth=month
gen double quarter=cond(missing(day),day,qofd(day)-_n)
gen double leftquarter=quarter
gen double week=cond(missing(day),day,wofd(day)+13*_n)
gen double leftweek=week
gen double half=cond(missing(day),day,hofd(day)+_n)
gen double lefthalf=half
format month %tm
format leftmonth %-tm
format quarter %tq
format leftquarter %-tqCCYY!qq
format week %tw
format leftweek %-tw
format half %th
format lefthalf %-th
foreach pair in "spss sav" "sasxport8 v8xpt" {
    gettoken kind ext : pair
    local ext=strtrim("`ext'")
    cio_export_test, kind(`kind') ext(`ext') name("`kind' period dates use first-day date formats") select(month leftmonth quarter leftquarter week leftweek half lefthalf)
}
cio_export_test, kind(sasxport8) ext(v8xpt) name("sasxport8 leap-second clocks") select(leap leftleap)
restore
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
cio_import_test using "temp/io_allformats/nonint.sav", kind(spss) name("SPSS non-integer and out-of-range value labels are skipped")

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
/* Text from a file (CSV/TSV headers, XLSX sheet names and header cells, XPORT
   labels, .dta variable labels exported as headers) is data. Each payload
   would write fio_m1.txt or fio_m2.txt if it ran as Stata code; one closes
   compound quotes, one is a macro expression, one holds $ and local
   references. The ctools command must give native's rc, variable names,
   labels, values and r() text, and create no marker. The payloads are built
   and compared in Mata only, so the test itself never expands them. */
capture mata: mata drop _fio_meta()
mata:
string rowvector _fio_meta(real scalar names)
{
    real scalar j
    string rowvector m
    m = J(1, 0, "")
    for (j = 1; j <= st_nvar(); j++) {
        m = m, (names ? st_varname(j) : ""), st_varlabel(j),
            (st_isstrvar(j) ? invtokens(st_sdata(., j)') : invtokens(strofreal(st_data(., j))'))
    }
    return(m)
}
_fio_q = char(34); _fio_b = char(96); _fio_a = char(39)
_fio_p = ("x" + _fio_q + _fio_a + "+string(filewrite(" + _fio_q + "fio_m1.txt" + _fio_q + "," + _fio_q + "x" + _fio_q + "))+" + _fio_b + _fio_q + "y",
    _fio_b + "=filewrite(" + _fio_q + "fio_m2.txt" + _fio_q + "," + _fio_q + "x" + _fio_q + ")" + _fio_a,
    "cost " + char(36) + "HOME " + _fio_b + "c(os)" + _fio_a + " " + _fio_q + "q" + _fio_q)
end
capture program drop _fio_case
program define _fio_case
    syntax, native(string asis) ctools(string asis) name(string) [rmacros(string) NONAMES]
    foreach side in native ctools {
        clear
        capture ``side''
        local rc_`side' = _rc
        mata: _fio_r_`side' = J(1, 0, "")
        foreach r of local rmacros {
            mata: _fio_r_`side' = _fio_r_`side', st_global("r(`r')")
        }
        mata: _fio_m_`side' = _fio_meta("`nonames'" == "")
    }
    local rc = 0
    if `rc_native' != `rc_ctools' local rc = 9
    mata: st_local("same", strofreal(_fio_m_native == _fio_m_ctools & _fio_r_native == _fio_r_ctools))
    if !`same' local rc = 9
    foreach marker in fio_m1.txt fio_m2.txt {
        capture confirm file "`marker'"
        if !_rc local rc = 9
        capture erase "`marker'"
    }
    if `rc' di as error "`name': rc native `rc_native' ctools `rc_ctools', same `same'"
    cio_record "`name'" `rc'
end
* A dataset whose variable labels carry the payloads, exported with them as
* the header row, then read back as text by native import excel
capture program drop _fio_export
program define _fio_export
    args command file
    clear
    quietly set obs 2
    foreach v in a b c {
        quietly gen double `v' = _n
    }
    mata: st_varlabel("a", _fio_p[1]); st_varlabel("b", _fio_p[2]); st_varlabel("c", _fio_p[3])
    `command' excel using "`file'", firstrow(varlabels) replace
    import excel using "`file'", allstring clear
end
_fio_case, native(import delimited using "temp/io_allformats/inj_header.csv", clear) ///
    ctools(cimport delimited using "temp/io_allformats/inj_header.csv", clear) ///
    name("file text: CSV header cells")
_fio_case, native(import delimited using "temp/io_allformats/inj_header.tsv", clear delimiters(tab)) ///
    ctools(cimport delimited using "temp/io_allformats/inj_header.tsv", clear delimiters(tab)) ///
    name("file text: TSV header cell")
_fio_case, native(import excel using "temp/io_allformats/inj_text.xlsx", describe) ///
    ctools(cimport excel using "temp/io_allformats/inj_text.xlsx", describe) ///
    name("file text: XLSX sheet name in describe") rmacros(worksheet_1 range_1)
_fio_case, native(import excel using "temp/io_allformats/inj_text.xlsx", firstrow clear) ///
    ctools(cimport excel using "temp/io_allformats/inj_text.xlsx", firstrow clear) ///
    name("file text: XLSX header cells")
clear
quietly set obs 2
gen double a = _n
gen double b = _n
mata: st_varlabel("a", _fio_p[2]); st_varlabel("b", _fio_p[3])
export sasxport5 a b using "temp/io_allformats/inj_label.xpt", replace
_fio_case, native(import sasxport5 "temp/io_allformats/inj_label.xpt", describe) ///
    ctools(cimport sasxport5 "temp/io_allformats/inj_label.xpt", describe) ///
    name("file text: XPORT labels in describe") rmacros(members)
_fio_case, native(import sasxport5 "temp/io_allformats/inj_label.xpt", clear) ///
    ctools(cimport sasxport5 "temp/io_allformats/inj_label.xpt", clear) ///
    name("file text: XPORT labels")
_fio_case, native(_fio_export export "temp/io_allformats/native/inj_header.xlsx") ///
    ctools(_fio_export cexport "temp/io_allformats/result/inj_header.xlsx") ///
    name("file text: variable labels as exported headers")
mata: mata drop _fio_q _fio_b _fio_a _fio_p
capture mata: mata drop _fio_r_native _fio_r_ctools _fio_m_native _fio_m_ctools

* SPSS: a string variable's value labels are skipped (native notes it), and
* label sets are named labels# by their dictionary ordinal, not by use.
cio_import_test using "temp/io_allformats/strlab.sav", kind(spss) name("SPSS string variable with value labels")
cio_import_test using "temp/io_allformats/order.sav", kind(spss) name("SPSS label set names follow the dictionary")
* XPORT 5: value labels are named after the format name alone, never after
* the SAS formats native knows; a string variable with a format fails like
* native unless novallabels is given.
foreach f in fmtnum fmtuser fmtnames {
    foreach opt in "" "novallabels" {
        cio_import_test using "temp/io_allformats/`f'.xpt", kind(sasxport5) name("XPORT5 format names as value labels `f' `opt'") opts(`opt')
    }
}
cio_import_test using "temp/io_allformats/fmtchar.xpt", kind(sasxport5) name("XPORT5 string format with novallabels") opts(novallabels)
capture noisily import sasxport5 "temp/io_allformats/fmtchar.xpt", clear
local nr=_rc
capture noisily cimport sasxport5 "temp/io_allformats/fmtchar.xpt", clear
cio_record "XPORT5 string format without novallabels fails as native" `=cond(`nr'==610 & _rc==610,0,9)'
capture mkdir "temp/io_allformats/xpf"
copy "temp/io_allformats/fmtuser.xpt" "temp/io_allformats/xpf/fmtuser.xpt", replace
clear
quietly set obs 2
gen yesno = _n-1
gen foo = _n
label define yesno 0 "No" 1 "Yes"
label define foo 1 "one" 2 "two"
label values yesno yesno
label values foo foo
export sasxport5 "temp/io_allformats/xpf/labels.xpt", replace vallabfile(xpf)
cio_import_test using "temp/io_allformats/xpf/fmtuser.xpt", kind(sasxport5) name("XPORT5 YESNO5. takes the YESNO labels of formats.xpf")
* dBase numeric fields holding non-integers import as double, as native,
* even when float holds each value exactly.
clear
quietly set obs 3
gen float quarter = _n/4
gen double half = _n + .5
gen float tenth = 1.1*_n
export dbase using "temp/io_allformats/native/exactfloat.dbf", replace
cio_import_test using "temp/io_allformats/native/exactfloat.dbf", kind(dbase) name("dBase exact-float non-integers are double")

/* XPORT 5 value-label names that are not SAS format names as written: two
   characters, trailing digits (native writes LAB1 as LAB). */
foreach nm in yn dt _x x_ gl ab1 ab12 lab1 a2 x9y9 {
    clear
    quietly set obs 2
    gen x = _n
    label define `nm' 1 "one" 2 "two"
    label values x `nm'
    cio_export_test, kind(sasxport5) ext(xpt) name("XPORT5 value label named `nm'")
}
/* File names that are not SAS names: native names a version 5 member after
   the file (to its first period, eight bytes, as written) and a version 8
   member DATASET. */
capture mkdir "temp/io_allformats/names"
capture mkdir "temp/io_allformats/names/native"
capture mkdir "temp/io_allformats/names/result"
clear
quietly set obs 2
gen x = _n
foreach b in "a-b" "1a" "a.b" "x y" "results-2020" "_N_" "a+b" "abcdefghijk" {
    foreach pair in "sasxport5 xpt" "sasxport8 v8xpt" {
        gettoken kind ext : pair
        local ext = strtrim("`ext'")
        local clause = cond("`kind'"=="sasxport5", "", "using")
        capture noisily export `kind' x using "temp/io_allformats/names/native/`b'.`ext'", replace
        local nr = _rc
        capture noisily cexport `kind' x using "temp/io_allformats/names/result/`b'.`ext'", replace
        local cr = _rc
        local rc = cond(`nr' == `cr', 0, 9)
        if !`rc' & !`nr' {
            preserve
            quietly import `kind' `clause' "temp/io_allformats/names/native/`b'.`ext'", clear
            tempfile nimp
            quietly save `nimp'
            quietly import `kind' `clause' "temp/io_allformats/names/result/`b'.`ext'", clear
            capture cf _all using `nimp', all
            if _rc local rc = 9
            restore
            if "`kind'" == "sasxport5" {
                quietly import sasxport5 "temp/io_allformats/names/native/`b'.`ext'", describe
                local nm `"`r(members)'"'
                quietly import sasxport5 "temp/io_allformats/names/result/`b'.`ext'", describe
                if `"`r(members)'"' != `"`nm'"' local rc = 9
            }
        }
        cio_record "`kind' file name `b' (native rc=`nr', C rc=`cr')" `rc'
    }
}
/* A write-protected destination is not replaced: native's code per writer,
   and the file is unchanged. Native SPSS and XPORT 8 delete their own file
   before failing; cexport keeps it. */
local ro = cond(c(os) == "Windows", "attrib +r", "chmod 444")
local rw = cond(c(os) == "Windows", "attrib -r", "chmod 644")
sysuse auto, clear
foreach pair in "delimited csv" "excel xlsx" "excel xls" "spss sav" "sasxport5 xpt" "sasxport8 v8xpt" "dbase dbf" {
    gettoken kind ext : pair
    local ext = strtrim("`ext'")
    foreach side in native result {
        local cmd = cond("`side'" == "native", "export", "cexport")
        local f "temp/io_allformats/`side'/readonly.`ext'"
        capture shell `rw' "`f'"
        quietly `cmd' `kind' make price using "`f'", replace
        shell `ro' "`f'"
        quietly checksum "`f'"
        local before = r(checksum)
        capture noisily `cmd' `kind' make using "`f'", replace
        local rc_`side' = _rc
        capture checksum "`f'"
        local same_`side' = !_rc & r(checksum) == `before'
        capture shell `rw' "`f'"
    }
    cio_record "`kind' `ext' read-only destination (native rc=`rc_native', C rc=`rc_result')" `=cond(`rc_native' == `rc_result' & `rc_native' & `same_result', 0, 9)'
}
import shp using "temp/io_allformats/shape1.shp", clear
foreach side in native result {
    local cmd = cond("`side'" == "native", "export", "cexport")
    local f "temp/io_allformats/`side'/readonly.shp"
    capture shell `rw' "`f'"
    quietly `cmd' shp using "`f'", replace
    shell `ro' "`f'"
    quietly checksum "`f'"
    local before = r(checksum)
    capture noisily `cmd' shp using "`f'", replace
    local rc_`side' = _rc
    capture checksum "`f'"
    local same_`side' = !_rc & r(checksum) == `before'
    capture shell `rw' "`f'"
}
cio_record "shp read-only destination (native rc=`rc_native', C rc=`rc_result')" `=cond(`rc_native' == `rc_result' & `rc_native' & `same_result', 0, 9)'
/* replace onto a symbolic link replaces the link with a regular file and
   leaves the target; sheet(..., modify) writes through the link. */
if c(os) != "Windows" {
    foreach pair in "delimited csv" "excel xlsx" "excel xls" {
        sysuse auto, clear
        gettoken kind ext : pair
        local ext = strtrim("`ext'")
        foreach side in native result {
            local cmd = cond("`side'" == "native", "export", "cexport")
            local d "temp/io_allformats/`side'"
            capture erase "`d'/link.`ext'"
            quietly `cmd' `kind' make using "`d'/target.`ext'", replace
            quietly checksum "`d'/target.`ext'"
            local target_`side' = r(checksum)
            shell ln -s target.`ext' "`d'/link.`ext'"
            capture noisily `cmd' `kind' price using "`d'/link.`ext'", replace
            local rc_`side' = _rc
            quietly checksum "`d'/target.`ext'"
            local kept_`side' = r(checksum) == `target_`side''
        }
        local rc = cond(`rc_native' == `rc_result' & `kept_native' & `kept_result', 0, 9)
        if !`rc' {
            if "`kind'" == "delimited" {
                quietly import delimited using "temp/io_allformats/native/link.`ext'", clear
                tempfile nimp
                quietly save `nimp'
                quietly import delimited using "temp/io_allformats/result/link.`ext'", clear
            }
            else {
                quietly import excel using "temp/io_allformats/native/link.`ext'", clear
                tempfile nimp
                quietly save `nimp'
                quietly import excel using "temp/io_allformats/result/link.`ext'", clear
            }
            capture cf _all using `nimp', all
            if _rc local rc = 9
        }
        cio_record "`kind' `ext' replace onto a symbolic link (native rc=`rc_native', C rc=`rc_result')" `rc'
        if "`kind'" == "excel" {
            sysuse auto, clear
            foreach side in native result {
                local cmd = cond("`side'" == "native", "export", "cexport")
                local d "temp/io_allformats/`side'"
                capture erase "`d'/mlink.`ext'"
                quietly `cmd' `kind' make using "`d'/mtarget.`ext'", replace
                shell ln -s mtarget.`ext' "`d'/mlink.`ext'"
                capture noisily `cmd' `kind' price using "`d'/mlink.`ext'", sheet("Sheet1", modify) cell(C3)
                local rc_`side' = _rc
            }
            local rc = cond(`rc_native' == `rc_result', 0, 9)
            if !`rc' {
                quietly import excel using "temp/io_allformats/native/mtarget.`ext'", clear
                tempfile nimp
                quietly save `nimp'
                quietly import excel using "temp/io_allformats/result/mtarget.`ext'", clear
                capture cf _all using `nimp', all
                if _rc local rc = 9
                quietly import excel using "temp/io_allformats/result/mlink.`ext'", clear
                capture cf _all using `nimp', all
                if _rc local rc = 9
            }
            cio_record "`ext' sheet modify through a symbolic link (native rc=`rc_native', C rc=`rc_result')" `rc'
        }
    }
}
/* dBase: the record of 200 doubles is 4001 bytes; native limits c(width)
   (1600 here) to 4000, and limits names, strings and variable counts. */
clear
quietly set obs 3
forvalues i = 1/200 {
    quietly gen double v`i' = _n + `i'
}
cio_export_test, kind(dbase) ext(dbf) name("dBase 200 doubles (record 4001 bytes)")
foreach limit in width count {
    clear
    quietly set obs 2
    if "`limit'" == "width" {
        forvalues i = 1/16 {
            quietly gen str254 s`i' = "x"
        }
    }
    else {
        forvalues i = 1/256 {
            quietly gen byte v`i' = 1
        }
    }
    capture noisily export dbase using "temp/io_allformats/native/limit.dbf", replace
    local nr = _rc
    capture noisily cexport dbase using "temp/io_allformats/result/limit.dbf", replace
    cio_record "dBase `limit' limit (native rc=`nr', C rc=`=_rc')" `=cond(`nr' == _rc & `nr', 0, 9)'
}
/* SPSS display formats: width plus three, COMMA for comma formats, E for %e,
   decimals after a period or a decimal comma. */
clear
quietly set obs 2
local i 0
foreach f in %1.0f %5.1f %9.2f %12.2fc %12.0gc %10.3e %8.2gc %-9.2f %09.2f %9,2f %9,2fc %10.4g %20.0fc %7.0g %8.0g {
    local ++i
    quietly gen double v`i' = 1234.5678 * _n
    format v`i' `f'
}
cio_export_test, kind(spss) ext(sav) name("SPSS numeric display formats")

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
    syntax, NAME(string) [OPTS(string asis) SELECT(string asis)]
    quietly export delimited `select' using "temp/io_parity/reference.csv", replace `opts'
    capture noisily cexport delimited `select' using "temp/io_parity/result.csv", replace `opts'
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
* Compare bytes: auto-detecting an import can conceal the wrong delimiter.
local tab=char(9)
foreach delimiter in " " "`tab'" "tab" "," ";" {
    io_export_case, name("literal delimiter bytes") opts(delimiter("`delimiter'"))
    io_export_case, name("literal delimiter quoted selected bytes") opts(delimiter("`delimiter'") novarnames quote) select(d s if _n>1 in 1/5)
    quietly export delimited d s using "temp/io_parity/reference.csv", delimiter("`delimiter'") replace
    capture erase "temp/io_parity/new.csv"
    cexport delimited d s using "temp/io_parity/new.csv", delimiter("`delimiter'") threads(2)
    mata: assert(io_samefile("temp/io_parity/reference.csv","temp/io_parity/new.csv"))
}
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
